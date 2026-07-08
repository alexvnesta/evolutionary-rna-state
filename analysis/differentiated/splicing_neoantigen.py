"""
analysis/differentiated/splicing_neoantigen.py

DIFFERENTIATED antigen module — splicing-derived neoantigen burden.

Reference tool
--------------
SNAF — Splicing Neo Antigen Finder (Li et al. 2024, Sci. Transl. Med.,
doi:10.1126/scitranslmed.ade2886; github.com/frankligy/SNAF). SNAF's pipeline
is: per-sample splice-junction counts -> tumor-specific "neojunctions"
(junctions abundant in tumor but ~absent in a normal-tissue reference) ->
3-frame in-silico translation ACROSS the junction -> junction-spanning k-mer
peptides -> MHC-I binding prediction.

Install status on this box (arm64 macOS, Python 3.11)
-----------------------------------------------------
`pip install snaf` is NOT installable here: SNAF 0.7.0 (and every published
version 0.5.0-0.7.0) hard-pins ``tensorflow==2.3.0``, for which no
osx-arm64 / py3.11 wheel exists (`pip install --dry-run snaf` ->
ResolutionImpossible). SNAF is also Linux-oriented and expects an AltAnalyze
exon-junction database + a GTEx junction control h5ad that are not present in
the sandbox.

Per the task's documented fallback, this module WRAPS THE SNAF ALGORITHM
FAITHFULLY rather than importing the unbuildable package. Every algorithmic
choice below is a direct port of SNAF's source (snaf/snaf.py, snaf/gtex.py,
v0.7.0), with the divergences noted explicitly:

  SNAF step                     | this module
  ------------------------------|-----------------------------------------------
  crude_tumor_specificity()     | call_neojunctions(): keep junction iff
   gtex.py:252                  |   count - normal_mean >= t_min (20) AND
                                |   normal_mean < n_max (3). Identical logic.
  NeoJunction.in_silico_        | translate_junction(): 3-frame (phase 0/1/2)
   translation() snaf.py:1033   |   read-through of first exon into second exon.
  get_peptides() snaf.py:1186   | _get_peptides(): verbatim port — translate
                                |   first-to-stop, continue into second, emit
                                |   only k-mers that SPAN the junction.
  subexon_tran()/query_from_    | retrieve_flanking_seq(): SNAF reads flanking
   dict_fa() snaf.py:1308/1228  |   EXON sequence from AltAnalyze's DB; we read
                                |   it from the GRCh38 FASTA + GENCODE exon
                                |   annotation instead (no AltAnalyze DB here).
                                |   SUBSTITUTION — see note in the function.
  run_MHCflurry() binding.py    | analysis.antigen_core.mhc_binding.count_binders
                                |   (SHARED engine; MHCflurry 2.x, rank<=2.0).

So the peptide set is derived by SNAF's exact translation logic; only the
sequence-lookup backend and the binding backend are swapped for this project's
shared, reproducible infrastructure.

Input (contract)
----------------
Per-sample splice-junction COUNTS. SNAF's neojunction gate is count-based, so
the natural input is STAR's ``SJ.out.tab`` (emitted per sample by the pipeline
session's STAR rnaseq spine): the universal, per-sample junction-count file.
A generic tidy junction table (from rMATS ``fromGTF.*.txt`` / SUPPA2 ``.ioe``
coordinates + counts) is also accepted via ``junctions_from_frame``.

Output (contract)
-----------------
``splice_neoantigen_burden`` — int, per (run_accession, cohort). Unique MHC-I
binding peptides (rank<=2.0) among peptides translated from tumor-specific
neojunctions. Comparable to every sibling ``*_neoantigen_burden`` because it
shares one binding engine.

Batch robustness
----------------
The neojunction gate uses raw read counts, which ARE library-size sensitive;
we therefore (a) express the tumor cutoff as a count threshold that the pilot
sets per platform, and (b) recommend a within-sample-normalized count (CPM) or
a shared GTEx-style normal reference so the gate is comparable across depth.
The binder step is the batch-invariant part: MHCflurry percentile rank is
allele-calibrated against a fixed random-peptide background, so it does not
drift with input composition. Report burden per clinical context, not pooled
across cohorts. (Same discipline as the sibling burden features.)
"""
from __future__ import annotations

import bisect
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
from Bio.Seq import Seq

from analysis.antigen_core.mhc_binding import (
    count_binders,
    clean_peptides,
    MIN_PEPTIDE_LEN,
    MAX_PEPTIDE_LEN,
)

# ---------------------------------------------------------------------------
# Defaults — ported from SNAF's initialize() (snaf/__init__.py)
# ---------------------------------------------------------------------------
T_MIN = 20      # SNAF t_min: min (count - normal_mean) for a tumor junction
N_MAX = 3       # SNAF n_max: max normal-tissue mean for a "neo" junction
# SNAF default in_silico_translation ks=[9,10]; we widen to the class-I window
# the shared engine scores (8-11) so the burden uses the full 8-11mer set.
KS = (8, 9, 10, 11)
DEFAULT_FLANK = 60      # nt of flanking exon sequence to read on each side.
                        # 60 nt -> 20 codons, enough to place any 8-11mer that
                        # spans the junction while keeping translation cheap.

_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FASTA = _REPO_ROOT / "reference" / "GRCh38" / "GRCh38.primary_assembly.genome.fa"
DEFAULT_GTF = _REPO_ROOT / "reference" / "GRCh38" / "gencode.v46.primary_assembly.annotation.gtf"


# ---------------------------------------------------------------------------
# Junction representation
# ---------------------------------------------------------------------------
@dataclass
class Junction:
    """A single splice junction with its per-sample read count.

    Coordinates follow the STAR ``SJ.out.tab`` convention: ``start``/``end`` are
    the 1-based first/last base of the INTRON (i.e. the junction gap). The
    donor exon ends at ``start-1``; the acceptor exon starts at ``end+1``.
    """
    chrom: str
    start: int          # 1-based first intronic base
    end: int            # 1-based last intronic base
    strand: str         # '+' or '-'
    count: int
    annotated: int = 0  # STAR col 6: 0 = novel, 1 = in annotation
    normal_mean: float = 0.0  # mean count in a normal-tissue reference (0 if none)

    @property
    def uid(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}({self.strand})"


# ---------------------------------------------------------------------------
# Input adapters
# ---------------------------------------------------------------------------
_STRAND_MAP = {"0": ".", "1": "+", "2": "-", 0: ".", 1: "+", 2: "-"}


def parse_star_sj(path: str | Path,
                  min_reads: int = 1,
                  novel_only: bool = False) -> list[Junction]:
    """Parse a STAR ``SJ.out.tab`` into a list of Junctions (SNAF-native input).

    STAR ``SJ.out.tab`` columns (tab-separated, no header):
        1 chrom  2 intron_start(1-based)  3 intron_end(1-based)
        4 strand (0 undef, 1 +, 2 -)      5 intron motif   6 annotated (0/1)
        7 n_uniquely_mapping_reads        8 n_multimapping_reads  9 max_overhang

    The junction COUNT used is the uniquely-mapping read count (col 7) — the
    quantity SNAF's tumor-specificity gate consumes.

    Parameters
    ----------
    min_reads : drop junctions with fewer than this many unique reads (noise).
    novel_only : keep only junctions STAR flags as not-in-annotation (col 6 == 0).
        A cheap first-pass proxy for "candidate neojunction" when no normal
        reference is supplied; the count gate in ``call_neojunctions`` is the
        real filter.
    """
    junctions: list[Junction] = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            strand = _STRAND_MAP.get(f[3], ".")
            if strand == ".":
                continue  # strand-undefined junctions are not translatable
            n_unique = int(f[6])
            if n_unique < min_reads:
                continue
            annotated = int(f[5]) if len(f) > 5 else 0
            if novel_only and annotated != 0:
                continue
            junctions.append(Junction(
                chrom=f[0], start=int(f[1]), end=int(f[2]),
                strand=strand, count=n_unique, annotated=annotated,
            ))
    return junctions


def junctions_from_frame(df: pd.DataFrame,
                         count_col: str = "count",
                         normal_mean_col: str | None = None) -> list[Junction]:
    """Adapter for a generic tidy junction table (rMATS/SUPPA2-derived).

    Expects columns ``chrom, start, end, strand, <count_col>`` where start/end
    are the 1-based intron boundaries. ``normal_mean_col`` (optional) supplies a
    per-junction normal-tissue mean count for the SNAF gate.
    """
    out: list[Junction] = []
    for _, r in df.iterrows():
        out.append(Junction(
            chrom=str(r["chrom"]), start=int(r["start"]), end=int(r["end"]),
            strand=str(r["strand"]), count=int(r[count_col]),
            annotated=int(r.get("annotated", 0)) if "annotated" in df.columns else 0,
            normal_mean=float(r[normal_mean_col]) if normal_mean_col else 0.0,
        ))
    return out


# ---------------------------------------------------------------------------
# Step 1 — neojunction calling  (port of gtex.crude_tumor_specificity)
# ---------------------------------------------------------------------------
def is_neojunction(count: float, normal_mean: float,
                   t_min: int = T_MIN, n_max: int = N_MAX) -> bool:
    """SNAF crude tumor-specificity test (gtex.py:252), verbatim logic.

    A junction is a neojunction iff it is abundant in tumor beyond the normal
    mean by at least ``t_min`` AND the normal mean itself is below ``n_max``::

        (count - normal_mean) >= t_min  AND  normal_mean < n_max

    With no normal reference (``normal_mean == 0``) this reduces to
    ``count >= t_min`` — the correct degenerate behaviour.
    """
    return (normal_mean < n_max) and ((count - normal_mean) >= t_min)


def call_neojunctions(junctions: Sequence[Junction],
                      t_min: int = T_MIN, n_max: int = N_MAX) -> list[Junction]:
    """Filter junctions to tumor-specific neojunctions (SNAF gate)."""
    return [j for j in junctions if is_neojunction(j.count, j.normal_mean, t_min, n_max)]


# ---------------------------------------------------------------------------
# Exon annotation index (replaces SNAF's AltAnalyze exon DB)
# ---------------------------------------------------------------------------
@dataclass
class ExonIndex:
    """Per-(chrom, strand) sorted exon-boundary index built from a GTF.

    SNAF looks flanking exon sequence up in AltAnalyze's exon-junction DB. We
    don't have that DB, so we index GENCODE exon boundaries: for the donor side
    we need exons ENDING at ``intron_start-1``; for the acceptor side, exons
    STARTING at ``intron_end+1``. Sequence itself comes from the FASTA.
    """
    # (chrom, strand) -> {exon_end -> [exon_start,...]}  (donor lookup)
    ends: dict = field(default_factory=dict)
    # (chrom, strand) -> {exon_start -> [exon_end,...]}  (acceptor lookup)
    starts: dict = field(default_factory=dict)

    @classmethod
    def from_gtf(cls, gtf_path: str | Path, chroms: set[str] | None = None) -> "ExonIndex":
        """Build the index from a GTF, optionally restricted to ``chroms``.

        Restricting to the chromosomes actually present in the junction set
        keeps memory/time bounded (the full GENCODE GTF is ~1.6 GB).
        """
        ends: dict = {}
        starts: dict = {}
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 8 or f[2] != "exon":
                    continue
                chrom = f[0]
                if chroms is not None and chrom not in chroms:
                    continue
                start, end, strand = int(f[3]), int(f[4]), f[6]
                key = (chrom, strand)
                ends.setdefault(key, {}).setdefault(end, []).append(start)
                starts.setdefault(key, {}).setdefault(start, []).append(end)
        return cls(ends=ends, starts=starts)


# ---------------------------------------------------------------------------
# Step 2 — flanking-sequence retrieval  (SUBSTITUTION for SNAF subexon_tran)
# ---------------------------------------------------------------------------
class FastaSeq:
    """Thin pysam FastaFile wrapper (builds the .fai index on first use)."""

    def __init__(self, fasta_path: str | Path):
        import pysam
        fasta_path = str(fasta_path)
        if not os.path.exists(fasta_path + ".fai"):
            pysam.faidx(fasta_path)
        self._fa = pysam.FastaFile(fasta_path)
        self.references = set(self._fa.references)

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        """0-based half-open fetch, returns uppercase sequence ('' if OOB)."""
        if chrom not in self.references:
            return ""
        start0 = max(0, start0)
        try:
            return self._fa.fetch(chrom, start0, end0).upper()
        except Exception:
            return ""


def retrieve_flanking_seq(j: Junction,
                          fasta: FastaSeq,
                          exon_index: ExonIndex | None = None,
                          flank: int = DEFAULT_FLANK) -> tuple[str, str]:
    """Return (first_seq, second_seq): donor- and acceptor-side exon sequence
    in TRANSCRIPTION orientation, ready for SNAF-style read-through translation.

    SUBSTITUTION NOTE
    -----------------
    SNAF pulls these two flanks from AltAnalyze's precomputed exon-sequence DB
    (``subexon_tran`` -> ``query_from_dict_fa``). We reconstruct the same two
    flanks from the GRCh38 FASTA:

      * If an ``ExonIndex`` is supplied and an annotated exon abuts the junction
        boundary, we take that exon's sequence (bounded to ``flank`` nt) — the
        faithful analogue of reading the annotated subexon.
      * Otherwise (novel boundary, or no index) we take ``flank`` nt of genomic
        sequence flanking the boundary. This matches SNAF's behaviour for novel
        exons / trailing coordinates, where it also reads raw flanking sequence.

    Orientation: for '+' strand, first = upstream (donor) flank ending at
    ``start-1``, second = downstream (acceptor) flank starting at ``end+1``.
    For '-' strand the roles swap and both flanks are reverse-complemented, so
    the returned pair is always 5'->3' in transcription order (as SNAF expects).
    """
    chrom, strand = j.chrom, j.strand
    donor_boundary = j.start - 1     # last exonic base before the intron (1-based)
    acceptor_boundary = j.end + 1    # first exonic base after the intron (1-based)

    def donor_flank_plus() -> str:
        # exon ending at donor_boundary; take its last `flank` bases
        if exon_index is not None:
            estarts = exon_index.ends.get((chrom, strand), {}).get(donor_boundary)
            if estarts:
                exon_start = max(estarts)  # shortest abutting exon -> tightest seq
                s = max(exon_start, donor_boundary - flank + 1)
                return fasta.fetch(chrom, s - 1, donor_boundary)
        return fasta.fetch(chrom, donor_boundary - flank, donor_boundary)

    def acceptor_flank_plus() -> str:
        if exon_index is not None:
            eends = exon_index.starts.get((chrom, strand), {}).get(acceptor_boundary)
            if eends:
                exon_end = min(eends)
                e = min(exon_end, acceptor_boundary + flank - 1)
                return fasta.fetch(chrom, acceptor_boundary - 1, e)
        return fasta.fetch(chrom, acceptor_boundary - 1, acceptor_boundary - 1 + flank)

    if strand == "+":
        first = donor_flank_plus()
        second = acceptor_flank_plus()
    else:
        # On '-' strand, transcription runs right->left. The donor (5') exon is
        # the one on the RIGHT (ending at acceptor_boundary side in genome), so
        # we fetch the genomic flanks and reverse-complement.
        # first (donor, transcription 5') = genomic sequence to the RIGHT of the
        #   intron (starting at end+1), reverse-complemented.
        if exon_index is not None:
            eends = exon_index.starts.get((chrom, strand), {}).get(acceptor_boundary)
            if eends:
                exon_end = min(eends)
                e = min(exon_end, acceptor_boundary + flank - 1)
                first_fwd = fasta.fetch(chrom, acceptor_boundary - 1, e)
            else:
                first_fwd = fasta.fetch(chrom, acceptor_boundary - 1, acceptor_boundary - 1 + flank)
            estarts = exon_index.ends.get((chrom, strand), {}).get(donor_boundary)
            if estarts:
                exon_start = max(estarts)
                s = max(exon_start, donor_boundary - flank + 1)
                second_fwd = fasta.fetch(chrom, s - 1, donor_boundary)
            else:
                second_fwd = fasta.fetch(chrom, donor_boundary - flank, donor_boundary)
        else:
            first_fwd = fasta.fetch(chrom, acceptor_boundary - 1, acceptor_boundary - 1 + flank)
            second_fwd = fasta.fetch(chrom, donor_boundary - flank, donor_boundary)
        first = str(Seq(first_fwd).reverse_complement()) if first_fwd else ""
        second = str(Seq(second_fwd).reverse_complement()) if second_fwd else ""
    return first, second


# ---------------------------------------------------------------------------
# Step 3 — in-silico translation  (verbatim port of SNAF get_peptides + loop)
# ---------------------------------------------------------------------------
def _get_peptides(de_facto_first: str, second: str, ks: Sequence[int],
                  phase: int) -> dict[int, list[str]]:
    """Verbatim port of SNAF snaf.get_peptides (snaf.py:1186), minus the
    start-codon-evidence bookkeeping we don't carry.

    Translate the first exon flank to a stop; if it reads through fully,
    continue the reading frame into the second exon flank and emit every
    k-mer that SPANS the junction (contains >=1 residue from the first part,
    except the in-frame extra==0 case which SNAF also requires >=1 from first).
    """
    out: dict[int, list[str]] = {k: [] for k in ks}
    extra = len(de_facto_first) % 3
    num = len(de_facto_first) // 3
    aa_first = str(Seq(de_facto_first).translate(to_stop=True))
    if len(aa_first) != num:          # premature stop in the first part -> abort
        return out
    if extra == 0:
        continue_second = second
    elif extra == 1:
        continue_second = de_facto_first[-1] + second
    else:  # extra == 2
        continue_second = de_facto_first[-2:] + second
    aa_second = str(Seq(continue_second).translate(to_stop=True))
    if len(aa_second) == 0:
        return out
    for k in ks:
        second_most = min(k, len(aa_second))
        first_most = len(aa_first)
        for n_from_second in range(second_most, 0, -1):
            n_from_first = k - n_from_second
            if n_from_first == 0 and extra == 0:
                # peptide entirely in second exon -> not splice-derived; skip
                continue
            if n_from_first <= first_most:
                if n_from_first > 0:
                    pep = aa_first[-n_from_first:] + aa_second[:n_from_second]
                else:
                    pep = aa_second[:n_from_second]
                out[k].append(pep)
    return out


def translate_junction(first: str, second: str,
                       ks: Sequence[int] = KS) -> list[str]:
    """SNAF 3-frame junction translation (NeoJunction.in_silico_translation).

    For each phase in {0,1,2}, drop the first ``phase`` bases of the donor flank
    and run ``_get_peptides``. Returns the deduped union of junction-spanning
    k-mer peptides across all three frames.
    """
    if not first or not second:
        return []
    peptides: set[str] = set()
    for phase in (0, 1, 2):
        de_facto_first = first[phase:]
        pep_dict = _get_peptides(de_facto_first, second, ks, phase)
        for k in ks:
            peptides.update(pep_dict[k])
    return sorted(peptides)


# ---------------------------------------------------------------------------
# Step 2+3 driver — neojunctions -> peptides
# ---------------------------------------------------------------------------
def peptides_from_neojunctions(neojunctions: Sequence[Junction],
                               fasta: FastaSeq,
                               exon_index: ExonIndex | None = None,
                               ks: Sequence[int] = KS,
                               flank: int = DEFAULT_FLANK) -> list[str]:
    """Translate a set of neojunctions into candidate 8-11mer peptides.

    Returns the deduped, hygiene-cleaned peptide list (standard-AA 8-11mers)
    ready for the shared binding engine.
    """
    peps: set[str] = set()
    for j in neojunctions:
        first, second = retrieve_flanking_seq(j, fasta, exon_index, flank=flank)
        peps.update(translate_junction(first, second, ks=ks))
    # engine hygiene (uppercases, dedupes, drops non-standard AA / out-of-range)
    return clean_peptides(peps)


# ---------------------------------------------------------------------------
# Top-level per-sample feature
# ---------------------------------------------------------------------------
def splice_neoantigen_burden(junctions: Sequence[Junction],
                             hla_alleles: Sequence[str],
                             fasta: FastaSeq | str | Path | None = None,
                             exon_index: ExonIndex | None = None,
                             t_min: int = T_MIN, n_max: int = N_MAX,
                             ks: Sequence[int] = KS,
                             flank: int = DEFAULT_FLANK,
                             rank_threshold: float = 2.0,
                             return_detail: bool = False):
    """Compute ``splice_neoantigen_burden`` for ONE sample.

    junction counts -> neojunctions (SNAF gate) -> junction-spanning peptides
    (SNAF 3-frame translation) -> SHARED count_binders -> burden int.

    Parameters
    ----------
    junctions : this sample's Junctions (from ``parse_star_sj`` /
        ``junctions_from_frame``), each carrying its read count (+ optional
        ``normal_mean``).
    hla_alleles : the sample's HLA-I alleles (from antigen_core.hla_typing).
    fasta : FastaSeq, or a path to the genome FASTA (defaults to the repo
        GRCh38). Required to translate — pass a small FASTA in tests.
    rank_threshold : binder percentile cutoff (default 2.0 = weak; pass 0.5 for
        strong-only), threaded straight through to the shared engine.

    Returns
    -------
    int  (the burden), or if ``return_detail`` a dict with the intermediate
    counts for QC/provenance.
    """
    if fasta is None:
        fasta = DEFAULT_FASTA
    if not isinstance(fasta, FastaSeq):
        fasta = FastaSeq(fasta)

    neojunctions = call_neojunctions(junctions, t_min=t_min, n_max=n_max)
    peptides = peptides_from_neojunctions(
        neojunctions, fasta, exon_index=exon_index, ks=ks, flank=flank)
    burden = count_binders(peptides, hla_alleles, rank_threshold=rank_threshold)

    if return_detail:
        return {
            "splice_neoantigen_burden": int(burden),
            "n_junctions": len(junctions),
            "n_neojunctions": len(neojunctions),
            "n_candidate_peptides": len(peptides),
        }
    return int(burden)


def build_feature_table(sample_junctions: dict[tuple[str, str], Sequence[Junction]],
                        sample_hla: dict[tuple[str, str], Sequence[str]],
                        fasta: FastaSeq | str | Path | None = None,
                        exon_index: ExonIndex | None = None,
                        **kwargs) -> pd.DataFrame:
    """Assemble the per-sample ``splice_neoantigen_burden`` feature table.

    Parameters
    ----------
    sample_junctions : {(run_accession, cohort): [Junction, ...]}
    sample_hla       : {(run_accession, cohort): [allele, ...]}

    Returns a tidy table keyed on (run_accession, cohort) per the v2 feature
    contract, with the single named feature column ``splice_neoantigen_burden``.
    """
    if fasta is None:
        fasta = DEFAULT_FASTA
    if not isinstance(fasta, FastaSeq):
        fasta = FastaSeq(fasta)
    rows = []
    for (run_accession, cohort), junctions in sample_junctions.items():
        alleles = sample_hla.get((run_accession, cohort), [])
        burden = splice_neoantigen_burden(
            junctions, alleles, fasta=fasta, exon_index=exon_index, **kwargs)
        rows.append({
            "run_accession": run_accession,
            "cohort": cohort,
            "splice_neoantigen_burden": int(burden),
        })
    return pd.DataFrame(rows, columns=["run_accession", "cohort",
                                       "splice_neoantigen_burden"])
