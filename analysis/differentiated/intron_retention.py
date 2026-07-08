"""
analysis/differentiated/intron_retention.py

DIFFERENTIATED BUCKET — Retained-intron load + cryptic-ORF neoantigens.

Two named, interpretable per-sample features on top of the pipeline session's
intron-retention quantification (IRFinder-S / featureCounts IR ratios):

    retained_intron_load    (int)  — how many introns are retained above a
                                      documented IR-ratio threshold in a sample
                                      (a "differentiated / RNA-processing
                                      dysregulation" phenotype: read-through of
                                      unspliced introns).
    ir_neoantigen_burden    (int)  — MHC-I binder count over CANDIDATE PEPTIDES
                                      translated from those retained introns
                                      (junction-spanning read-through peptides +
                                      purely intronic ORF peptides), scored by
                                      the SHARED MHCflurry engine so it is
                                      directly comparable to the splice / TE /
                                      fusion / SNV burdens.

Design (per FEATURE_CONTRACT_v2.md and feature_registry.json):
  * Upstream input: the pipeline's intron-retention output —
      - intron_retention.parquet  (tidy-wide: run_accession, cohort, <intron_id...> = IR ratio)
        or the per-sample long tables (run_accession, cohort, intron_id, IR_ratio)
      - introns.saf               (genomic coordinates of each intron_id, from
                                    pipelines/intron_retention/bin/make_intron_saf.py)
      - intron2gene.tsv           (intron_id -> gene_id, intron_length)
      - GRCh38 genome FASTA       (for translating retained introns)
      - the sample's HLA-I alleles (hla_typing.parquet / build_hla_table)
  * Shared engine: analysis.antigen_core.mhc_binding.count_binders — the SAME
    peptide scorer every antigen module uses. We NEVER touch MHCflurry directly.

Why "cryptic ORFs" from retained introns generate neoantigens
-------------------------------------------------------------
A retained intron is present in the mature mRNA. When translated, the ribosome
reads through the 5' splice site into intronic sequence in the host gene's
reading frame until it meets a premature termination codon. The residues after
the exon->intron junction are NOT in any canonical protein — they are a novel,
tumour-specific antigenic sequence. Two peptide classes arise:
  (1) JUNCTION-SPANNING read-through peptides — 8-11mers that straddle the
      exon/intron boundary (part canonical exon, part novel intronic frame).
  (2) INTRONIC ORF peptides — 8-11mers wholly inside an intronic ORF (an
      internal ATG..stop within the retained intron).
Both are tiled as 8-11mers and handed to the shared MHC-I engine.

Batch robustness
----------------
See BATCH_ROBUSTNESS_NOTE below and the module docstring of
`compute_retained_intron_load`. In short: the IR *ratio* is a within-sample
intron/exon density ratio (self-normalising for library size), but the *count*
of retained introns scales with sequencing depth (deeper libraries make more
introns evaluable and push more across the threshold). We therefore also emit
a depth-normalised fraction and a within-cohort z-score, cap the candidate
intron set for the burden, and rely on MHCflurry's allele-calibrated percentile
rank (batch-invariant) for the binder call.

NO real-cohort feature values are fabricated here. This module is a runnable,
unit-validated feature builder; the full-cohort run is deferred to the pilot.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

# ---------------------------------------------------------------------------
# Shared antigen core — import, do NOT reimplement.
# The module lives at analysis/antigen_core/mhc_binding.py. Support both an
# installed-package import (analysis.antigen_core...) and a direct-path import
# (when the antigen_core dir is on sys.path, as in its own unit tests).
# ---------------------------------------------------------------------------
try:  # package-style import (repo root on sys.path)
    from analysis.antigen_core.mhc_binding import (  # type: ignore
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
    from analysis.antigen_core.hla_typing import ALLELE_COLS  # type: ignore
except Exception:  # pragma: no cover - fallback for path-style import
    import sys
    _CORE = Path(__file__).resolve().parents[1] / "antigen_core"
    if str(_CORE) not in sys.path:
        sys.path.insert(0, str(_CORE))
    from mhc_binding import (  # type: ignore
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
    from hla_typing import ALLELE_COLS  # type: ignore

# ---------------------------------------------------------------------------
# Documented constants (all overridable at call time).
# ---------------------------------------------------------------------------
IR_RETAINED_THRESHOLD = 0.10   # IR ratio above which an intron is "retained".
#   Matches the pipeline summary default (compute_ir_ratio.py --high-ir-threshold
#   = 0.1) and the common IRFinder "IRratio" call for a retained intron.
EXON_FLANK_NT = 30             # 30 nt (10 codons) of upstream exon kept before
#   the 5' splice site, so junction-spanning 8-11mers can be tiled.
PEP_MIN, PEP_MAX = 8, 11       # MHC-I peptide length window.
MAX_INTRONS_FOR_BURDEN = 300   # cap on candidate introns per sample (highest IR
#   first) — bounds peptide count so a deep library cannot inflate the burden
#   purely through depth. Documented batch-robustness lever.
MAX_READTHROUGH_AA = 100       # stop walking a read-through ORF after this many
#   residues past the junction (guards against pathological long ORFs).

STOP = "*"
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Standard genetic code (DNA codon -> amino acid; '*' = stop).
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": STOP, "TAG": STOP,
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": STOP, "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

BATCH_ROBUSTNESS_NOTE = (
    "IR ratio = intron_read_density / (intron_read_density + exon_read_density) "
    "is a WITHIN-SAMPLE ratio, so each intron's IR value is self-normalising for "
    "library size and does not drift with total depth. The DERIVED per-sample "
    "COUNT of retained introns is, however, depth-sensitive: deeper libraries "
    "make more introns evaluable (host-gene exonic coverage passes the upstream "
    "min-gene-exon-count gate) and push more introns across the IR threshold. "
    "Mitigations built in here: (1) count only introns evaluable in the sample "
    "(non-NA IR, i.e. host gene had sufficient exonic coverage upstream); "
    "(2) emit retained_intron_fraction = load / n_introns_evaluated, which is "
    "depth-robust; (3) emit retained_intron_load_cohortz = within-cohort robust "
    "z-score of the fraction, so batches are never pooled on raw counts; "
    "(4) for ir_neoantigen_burden, cap candidate introns to the top "
    f"{MAX_INTRONS_FOR_BURDEN} by IR ratio per sample so a deep library cannot "
    "inflate the peptide pool by depth alone; (5) the binder call uses MHCflurry "
    "percentile rank, which is calibrated against a fixed per-allele background "
    "and is therefore batch/platform-invariant. Report per clinical context / "
    "z-scored within cohort; never pool raw counts across platforms."
)


# ===========================================================================
# Sequence utilities (stdlib translation — no BioPython dependency needed, and
# fully deterministic for unit tests).
# ===========================================================================
def revcomp(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def translate_frame(seq: str, frame: int) -> str:
    """Translate ``seq`` in reading ``frame`` (0/1/2). Unknown codons -> 'X'."""
    prot = []
    for i in range(frame, len(seq) - 2, 3):
        prot.append(_CODON_TABLE.get(seq[i:i + 3].upper(), "X"))
    return "".join(prot)


def _kmers(prot: str, kmin: int = PEP_MIN, kmax: int = PEP_MAX) -> list[str]:
    """All kmin..kmax substrings of a stop-free protein string."""
    out = []
    n = len(prot)
    for k in range(kmin, kmax + 1):
        for i in range(0, n - k + 1):
            out.append(prot[i:i + k])
    return out


def _kmers_spanning(prot: str, junction_aa: int,
                    kmin: int = PEP_MIN, kmax: int = PEP_MAX) -> list[str]:
    """kmin..kmax substrings of ``prot`` that STRADDLE the junction.

    ``junction_aa`` is the residue index in ``prot`` of the first amino acid
    encoded (wholly or partly) by intronic sequence. A peptide straddles the
    junction if it contains at least one residue before and at least one at/after
    ``junction_aa`` — i.e. start < junction_aa <= end.
    """
    out = []
    n = len(prot)
    for k in range(kmin, kmax + 1):
        # windows [i, i+k) with i < junction_aa <= i+k
        lo = max(0, junction_aa - k)
        hi = min(junction_aa - 1, n - k)
        for i in range(lo, hi + 1):
            out.append(prot[i:i + k])
    return out


# ===========================================================================
# Intron coordinate index (from the pipeline's introns.saf + genome FASTA).
# ===========================================================================
@dataclass
class IntronCoord:
    intron_id: str
    chrom: str
    strand: str
    # pure-intronic survivor sub-intervals, 1-based inclusive, genomic order:
    blocks: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.blocks)

    @property
    def end(self) -> int:
        return max(e for _, e in self.blocks)


def load_intron_saf(saf_path: str | Path) -> dict[str, IntronCoord]:
    """Parse introns.saf (GeneID/Chr/Start/End/Strand) into IntronCoord objects.

    A single intron_id (GeneID column) may span several rows (make_intron_saf
    fragments a pure-intronic region after exon masking); we collect them all.
    """
    coords: dict[str, IntronCoord] = {}
    df = pd.read_csv(saf_path, sep="\t")
    for _, r in df.iterrows():
        iid = str(r["GeneID"])
        c = coords.get(iid)
        if c is None:
            c = IntronCoord(iid, str(r["Chr"]), str(r["Strand"]))
            coords[iid] = c
        c.blocks.append((int(r["Start"]), int(r["End"])))
    for c in coords.values():
        c.blocks.sort()
    return coords


class GenomeFasta:
    """Thin wrapper over a genome FASTA giving 1-based inclusive slices.

    Uses pysam.FastaFile when available (indexed, memory-light — the pilot path
    on the full GRCh38), and falls back to an in-memory dict of sequences for
    tiny synthetic test FASTAs (no .fai needed).
    """

    def __init__(self, fasta_path: str | Path):
        self.path = str(fasta_path)
        self._fa = None
        self._mem: dict[str, str] | None = None
        try:
            import pysam  # type: ignore
            self._fa = pysam.FastaFile(self.path)
        except Exception:
            self._mem = _read_fasta_to_dict(self.path)

    def fetch(self, chrom: str, start1: int, end1: int) -> str:
        """1-based inclusive [start1, end1] on the FORWARD strand."""
        if start1 > end1:
            return ""
        if self._fa is not None:
            return self._fa.fetch(chrom, start1 - 1, end1).upper()
        assert self._mem is not None
        return self._mem.get(chrom, "")[start1 - 1:end1].upper()

    def close(self):
        if self._fa is not None:
            self._fa.close()


def _read_fasta_to_dict(path: str | Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# ===========================================================================
# Feature 1 — retained_intron_load
# ===========================================================================
def _ir_long_from_matrix(ir_matrix: pd.DataFrame) -> pd.DataFrame:
    """Melt a tidy-wide intron_retention matrix to long (run_accession, cohort,
    intron_id, IR_ratio)."""
    id_cols = [c for c in ("run_accession", "cohort") if c in ir_matrix.columns]
    feat_cols = [c for c in ir_matrix.columns if c not in id_cols]
    long = ir_matrix.melt(id_vars=id_cols, value_vars=feat_cols,
                          var_name="intron_id", value_name="IR_ratio")
    return long


def _robust_z(x: pd.Series) -> pd.Series:
    """Median/MAD robust z-score. Returns 0 where MAD==0 (or <2 samples)."""
    x = x.astype(float)
    med = x.median()
    mad = (x - med).abs().median()
    if not mad or pd.isna(mad):
        return pd.Series(0.0, index=x.index)
    return (x - med) / (1.4826 * mad)


def compute_retained_intron_load(
    ir_data: pd.DataFrame,
    threshold: float = IR_RETAINED_THRESHOLD,
    weighted: bool = True,
) -> pd.DataFrame:
    """Per-sample retained-intron load from IR ratios.

    Parameters
    ----------
    ir_data : either the tidy-wide intron_retention matrix (run_accession,
        cohort, <intron_id...> = IR ratio) OR a long table with columns
        (run_accession, cohort, intron_id, IR_ratio).
    threshold : IR ratio above which an intron is counted as retained.
    weighted : also compute a weighted sum (sum of IR ratios of retained
        introns) — a load that credits strongly-retained introns more.

    Returns
    -------
    DataFrame keyed on (run_accession, cohort) with columns:
        retained_intron_load        int    — # introns with IR >= threshold (NAMED FEATURE)
        n_introns_evaluated         int    — # introns with non-NA IR (denominator)
        retained_intron_fraction    float  — load / n_introns_evaluated (depth-robust)
        retained_intron_load_weighted float — sum of IR ratios of retained introns
        retained_intron_load_cohortz  float — within-cohort robust z of the fraction

    Batch note: retained_intron_load is the interpretable headline count;
    retained_intron_fraction and *_cohortz are the depth-/batch-robust forms to
    use when pooling cohorts (see BATCH_ROBUSTNESS_NOTE).
    """
    if "IR_ratio" in ir_data.columns and "intron_id" in ir_data.columns:
        long = ir_data[["run_accession", "cohort", "intron_id", "IR_ratio"]].copy()
    else:
        long = _ir_long_from_matrix(ir_data)

    long["IR_ratio"] = pd.to_numeric(long["IR_ratio"], errors="coerce")
    evaluable = long.dropna(subset=["IR_ratio"])

    def _agg(g: pd.DataFrame) -> pd.Series:
        retained = g["IR_ratio"] >= threshold
        n_eval = int(len(g))
        load = int(retained.sum())
        d = {
            "retained_intron_load": load,
            "n_introns_evaluated": n_eval,
            "retained_intron_fraction": (load / n_eval) if n_eval else float("nan"),
        }
        if weighted:
            d["retained_intron_load_weighted"] = float(g.loc[retained, "IR_ratio"].sum())
        return pd.Series(d)

    out = (evaluable.groupby(["run_accession", "cohort"], sort=False)
           .apply(_agg, include_groups=False).reset_index())

    # within-cohort robust z-score of the depth-normalised fraction
    out["retained_intron_load_cohortz"] = (
        out.groupby("cohort")["retained_intron_fraction"]
        .transform(_robust_z)
    )
    out["retained_intron_load"] = out["retained_intron_load"].astype(int)
    out["n_introns_evaluated"] = out["n_introns_evaluated"].astype(int)
    return out


# ===========================================================================
# Feature 2 — ir_neoantigen_burden  (cryptic-ORF peptides -> shared engine)
# ===========================================================================
def cryptic_orf_peptides(
    coord: IntronCoord,
    genome: GenomeFasta,
    exon_flank_nt: int = EXON_FLANK_NT,
    kmin: int = PEP_MIN,
    kmax: int = PEP_MAX,
    max_readthrough_aa: int = MAX_READTHROUGH_AA,
) -> list[str]:
    """Translate one retained intron into candidate cryptic-ORF peptides.

    Produces two peptide classes, tiled as kmin..kmax-mers:
      (A) JUNCTION-SPANNING read-through peptides — the ribosome reads the
          upstream exon in some frame and continues into the intron until the
          first stop; we keep the kmers that straddle the exon/intron boundary.
          The true exonic frame is unknown at this layer, so all 3 frames of the
          (exon_flank + intron) construct are translated (a superset; MHCflurry
          filtering downstream keeps only presented peptides).
      (B) INTRONIC ORF peptides — internal ATG..stop ORFs wholly inside the
          intron, in all 3 frames, tiled as kmers.

    Coordinates come from introns.saf; sequence from the genome FASTA. For a '-'
    strand intron the retrieved sequence is reverse-complemented so translation
    is 5'->3' on the coding strand, and the upstream exon flank is taken on the
    genomic 3' side (which is 5' in transcription).
    """
    chrom, strand = coord.chrom, coord.strand
    # intron sequence = concatenation of pure-intronic survivor blocks in
    # transcription order.
    if strand == "-":
        # transcription order is genomic 3'->5'; take blocks descending and revcomp
        blocks = sorted(coord.blocks, reverse=True)
        intron_seq = "".join(revcomp(genome.fetch(chrom, s, e)) for s, e in blocks)
        # upstream exon flank is immediately 3' of the intron in genomic coords
        gmax = max(e for _, e in coord.blocks)
        flank_fwd = genome.fetch(chrom, gmax + 1, gmax + exon_flank_nt)
        exon_flank = revcomp(flank_fwd)
    else:
        blocks = sorted(coord.blocks)
        intron_seq = "".join(genome.fetch(chrom, s, e) for s, e in blocks)
        gmin = min(s for s, _ in coord.blocks)
        exon_flank = genome.fetch(chrom, gmin - exon_flank_nt, gmin - 1)

    peptides: set[str] = set()
    if not intron_seq:
        return []

    # ---- (A) junction-spanning read-through, all 3 frames of exon+intron ----
    construct = exon_flank + intron_seq
    jn_nt = len(exon_flank)  # nt index of first intronic base in the construct
    for frame in (0, 1, 2):
        prot = translate_frame(construct, frame)
        if not prot:
            continue
        # residue index of the first codon that includes an intronic base:
        # codon r covers construct nt [frame + 3r, frame + 3r + 3); it is the
        # first to touch the intron when frame + 3r + 3 > jn_nt.
        junction_aa = max(0, -(-(jn_nt - frame) // 3))  # ceil((jn_nt-frame)/3)
        # walk from the junction to the first downstream stop (the read-through
        # ORF); keep the pre-stop segment around the junction.
        # Build the local protein window: from some residues before the junction
        # up to the first stop after it.
        # find first stop at/after junction_aa
        stop_idx = prot.find(STOP, junction_aa)
        seg_end = len(prot) if stop_idx == -1 else stop_idx
        seg_end = min(seg_end, junction_aa + max_readthrough_aa)
        # left context: back up to previous stop (keep in-frame exonic residues)
        prev_stop = prot.rfind(STOP, 0, junction_aa)
        seg_start = 0 if prev_stop == -1 else prev_stop + 1
        segment = prot[seg_start:seg_end]
        if len(segment) >= kmin:
            local_junction = junction_aa - seg_start
            peptides.update(_kmers_spanning(segment, local_junction, kmin, kmax))

    # ---- (B) intronic ORFs (ATG..stop) in all 3 frames ----
    for frame in (0, 1, 2):
        prot = translate_frame(intron_seq, frame)
        # split on stops; within each stop-free segment, take ORFs from each Met
        for seg in prot.split(STOP):
            m = seg.find("M")
            while m != -1:
                orf = seg[m:]
                if len(orf) >= kmin:
                    peptides.update(_kmers(orf, kmin, kmax))
                nxt = seg.find("M", m + 1)
                if nxt == -1 or nxt - m > 200:  # only a few Met starts per segment
                    break
                m = nxt

    return sorted(peptides)


def compute_ir_neoantigen_burden(
    ir_data: pd.DataFrame,
    intron_coords: dict[str, IntronCoord],
    genome: GenomeFasta,
    hla_by_sample: dict[str, Sequence[str]],
    threshold: float = IR_RETAINED_THRESHOLD,
    rank_threshold: float = WEAK_BINDER_RANK,
    max_introns: int = MAX_INTRONS_FOR_BURDEN,
) -> pd.DataFrame:
    """Per-sample intron-retention neoantigen burden via the SHARED engine.

    For each sample: take its retained introns (IR >= threshold, highest IR
    first, capped at ``max_introns``), translate each into cryptic-ORF candidate
    peptides, pool the unique peptides, and hand them + the sample's HLA-I
    alleles to ``count_binders`` (the one shared MHCflurry engine). The result
    is the count of UNIQUE binder peptides — directly comparable to
    splice_/te_/fusion_/snv_indel_neoantigen_burden.

    Parameters
    ----------
    hla_by_sample : {run_accession: [allele, ...]} (up to 6 HLA-I alleles).
        Samples absent from this map get NA burden (HLA not typed).
    rank_threshold : WEAK_BINDER_RANK (2.0) default; STRONG_BINDER_RANK (0.5)
        for strong-only.

    Returns
    -------
    DataFrame (run_accession, cohort, ir_neoantigen_burden, n_candidate_peptides,
    n_retained_introns_used). ir_neoantigen_burden is <NA> where HLA is missing.
    """
    if "IR_ratio" in ir_data.columns and "intron_id" in ir_data.columns:
        long = ir_data[["run_accession", "cohort", "intron_id", "IR_ratio"]].copy()
    else:
        long = _ir_long_from_matrix(ir_data)
    long["IR_ratio"] = pd.to_numeric(long["IR_ratio"], errors="coerce")

    rows = []
    for (samp, cohort), g in long.groupby(["run_accession", "cohort"], sort=False):
        retained = (g.dropna(subset=["IR_ratio"])
                    .query("IR_ratio >= @threshold")
                    .sort_values("IR_ratio", ascending=False)
                    .head(max_introns))
        alleles = hla_by_sample.get(samp)
        # pool candidate peptides across this sample's retained introns
        peps: set[str] = set()
        n_used = 0
        for iid in retained["intron_id"]:
            coord = intron_coords.get(str(iid))
            if coord is None:
                continue
            n_used += 1
            peps.update(cryptic_orf_peptides(coord, genome))
        pep_list = sorted(peps)
        if not alleles:
            burden = pd.NA
        else:
            burden = count_binders(pep_list, list(alleles),
                                   rank_threshold=rank_threshold)
        rows.append({
            "run_accession": samp,
            "cohort": cohort,
            "ir_neoantigen_burden": burden,
            "n_candidate_peptides": len(pep_list),
            "n_retained_introns_used": n_used,
        })
    out = pd.DataFrame(rows)
    if "ir_neoantigen_burden" in out:
        out["ir_neoantigen_burden"] = out["ir_neoantigen_burden"].astype("Int64")
    return out


# ===========================================================================
# Convenience: build both feature columns and write the contract matrix.
# ===========================================================================
def hla_map_from_table(hla_table: pd.DataFrame) -> dict[str, list[str]]:
    """{run_accession: [6 HLA-I alleles]} from an hla_typing.parquet-style table."""
    cols = [c for c in ALLELE_COLS if c in hla_table.columns]
    out: dict[str, list[str]] = {}
    for _, r in hla_table.iterrows():
        alleles = [str(r[c]) for c in cols if pd.notna(r[c]) and str(r[c]).strip()]
        out[str(r["run_accession"])] = alleles
    return out


def build_ir_features(
    ir_data: pd.DataFrame,
    intron_saf: str | Path,
    genome_fasta: str | Path,
    hla_table: pd.DataFrame,
    threshold: float = IR_RETAINED_THRESHOLD,
) -> pd.DataFrame:
    """End-to-end: retained_intron_load + ir_neoantigen_burden per sample.

    Returns a tidy feature frame keyed on (run_accession, cohort) carrying the
    two NAMED feature columns plus their companions/QC. This is what lands in
    results/features/ (merged into antigen_features.parquet per the v2 contract).
    """
    load = compute_retained_intron_load(ir_data, threshold=threshold)
    coords = load_intron_saf(intron_saf)
    genome = GenomeFasta(genome_fasta)
    try:
        hla_map = hla_map_from_table(hla_table)
        burden = compute_ir_neoantigen_burden(
            ir_data, coords, genome, hla_map, threshold=threshold)
    finally:
        genome.close()
    feat = load.merge(burden, on=["run_accession", "cohort"], how="outer")
    return feat


__all__ = [
    "IR_RETAINED_THRESHOLD", "EXON_FLANK_NT", "PEP_MIN", "PEP_MAX",
    "MAX_INTRONS_FOR_BURDEN", "BATCH_ROBUSTNESS_NOTE",
    "revcomp", "translate_frame",
    "IntronCoord", "load_intron_saf", "GenomeFasta",
    "compute_retained_intron_load",
    "cryptic_orf_peptides", "compute_ir_neoantigen_burden",
    "hla_map_from_table", "build_ir_features",
]
