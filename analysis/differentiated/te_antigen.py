"""
analysis/differentiated/te_antigen.py

DIFFERENTIATED BUCKET — Transposable-element / ERV-derived antigen burden.

Produces, per sample, the ``te_antigen_burden`` family of NAMED, INTERPRETABLE
features: the number of unique MHC-I binder peptides derived from the ORFs of
transcriptionally ACTIVE TE/ERV loci. Reactivation of endogenous
retroviruses (ERVs) and other LTR/LINE elements is a well-described source of
tumour-specific antigens invisible to a WES/annotation neoantigen pipeline —
the exact WES-blind RNA-state phenotype this project set out to test.

WHY THIS IS A REAL ANTIGEN SOURCE (not noise)
    - ERV/LTR elements retain gag/pol/env-like ORFs; when their LTR promoters
      are de-repressed in tumours they are transcribed and can be translated
      into peptides that reach MHC-I. HERV-E, HERV-K(HML-2) env/gag epitopes
      are documented T-cell targets in RCC and other tumours.
    - LINE-1 (L1) ORF1p/ORF2p are re-expressed in many carcinomas and are a
      further TE-derived antigen source; we translate any expressed locus's
      genomic sequence in all frames rather than assuming a canonical ORF.

PIPELINE INPUTS CONSUMED (v1 HANDOFF_CONTRACT / v2 FEATURE_CONTRACT)
    te_locus.parquet   — Telescope per-locus EM counts  (rows=run_accession,
                         cols=TE locus id, values=reassigned read count)
    te_family.parquet  — per-family aggregate counts (LINE/SINE/LTR/ERV/DNA)
    + a locus annotation (RepeatMasker/GENCODE-TE: locus_id -> repeat_class,
      repeat_family, chrom, start, end, strand) and the locus genomic
      sequences (extracted once from the genome FASTA for the expressed loci).

THE SHARED-ENGINE CONTRACT (do NOT reimplement peptide scoring)
    derive candidate 8-11mer peptides from the phenotype matrix
        -> get the sample's HLA-I alleles (analysis/antigen_core/hla_typing.py)
        -> analysis/antigen_core/mhc_binding.count_binders(peps, alleles)
        -> te_antigen_burden (int).
    All *_neoantigen_burden features share this one engine, so the TE burden is
    directly comparable to the splicing / IR / fusion / SNV burdens.

NAMED FEATURE OUTPUTS (interpretable — never an opaque embedding)
    te_antigen_burden            int   headline: unique MHC-I binder peptides
                                       from all active TE/ERV loci (rank<=2.0)
    te_antigen_burden_strong     int   same, strong threshold (rank<=0.5)
    te_antigen_burden_LINE       int   family-resolved (LINE-1 / L1 etc.)
    te_antigen_burden_SINE       int   family-resolved (Alu / MIR)
    te_antigen_burden_LTR        int   family-resolved (all LTR/ERV incl. ERV)
    te_antigen_burden_ERV        int   family-resolved (ERV subset of LTR:
                                       ERV1/ERVK/ERVL/ERVL-MaLR/HERV)
    te_antigen_n_expressed_loci  int   QC: loci passing the activity filter
    te_antigen_n_binder_loci     int   QC: loci contributing >=1 binder peptide
    te_antigen_top_locus         str   locus-level: locus with the most binders

BATCH ROBUSTNESS (project hard constraint) — see module-level note
    build_batch_robustness_note() and the FEATURE_CONTRACT entry. TE
    quantification is multimap-sensitive; the activity call is a WITHIN-sample
    CPM normalization, and the burden is an MHCflurry percentile-rank count
    (allele-calibrated, composition-invariant). REQUIRES that every sample's
    Telescope counts come from the SAME STAR multimap settings.

This module is a RUNNABLE, UNIT-VALIDATED module (see test_te_antigen.py) +
a Nextflow subworkflow stub (te_antigen.nf). It never fabricates real-cohort
per-sample feature values.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd

# --- import the shared antigen core (peptide-scoring engine) ---------------
# Works whether run from the repo root (as a package) or from a Nextflow work
# dir with analysis/ on sys.path. We import the module functions, never a
# private MHCflurry call — the engine is the single source of truth.
_ANALYSIS_DIR = Path(__file__).resolve().parents[1]        # .../analysis
if str(_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(_ANALYSIS_DIR))
try:  # package-style import (preferred)
    from analysis.antigen_core.mhc_binding import (  # type: ignore
        best_per_peptide, count_binders,
        STRONG_BINDER_RANK, WEAK_BINDER_RANK,
        MIN_PEPTIDE_LEN, MAX_PEPTIDE_LEN,
    )
except Exception:  # flat import (Nextflow work dir / tests)
    from antigen_core.mhc_binding import (  # type: ignore
        best_per_peptide, count_binders,
        STRONG_BINDER_RANK, WEAK_BINDER_RANK,
        MIN_PEPTIDE_LEN, MAX_PEPTIDE_LEN,
    )

# ---------------------------------------------------------------------------
# TE FAMILY TAXONOMY
# ---------------------------------------------------------------------------
# The five interpretable buckets the task asks for. RepeatMasker annotates each
# locus with a `repeat_class` (a.k.a. class/family, e.g. "LINE/L1",
# "SINE/Alu", "LTR/ERVK", "DNA/hAT-Charlie"). We collapse to a top-level FAMILY
# bucket and separately flag the ERV subset of LTR (the retrovirus-derived
# elements with coding gag/pol/env potential — the highest-value antigen
# source).
FAMILY_BUCKETS = ("LINE", "SINE", "LTR", "DNA", "OTHER")

# ERV = LTR loci whose repeat subfamily is a known endogenous-retrovirus group.
# Substring match against the (upper-cased) repeat family / class token.
_ERV_TOKENS = ("ERV", "HERV", "MALR", "MER", "HML", "MMTV", "GYPSY")


def classify_family(repeat_class: str) -> str:
    """Collapse a RepeatMasker `repeat_class` string to a top-level bucket.

    Accepts the compound "class/family" token RepeatMasker emits
    (e.g. 'LINE/L1', 'SINE/Alu', 'LTR/ERVK', 'DNA/TcMar-Tigger') or a bare
    class ('LINE'). Returns one of FAMILY_BUCKETS.
    """
    if repeat_class is None:
        return "OTHER"
    top = str(repeat_class).strip().upper().split("/", 1)[0]
    if top in ("LINE", "SINE", "LTR", "DNA"):
        return top
    # some annotations write ERV directly at top level -> treat as LTR bucket
    if top in ("ERV", "ERVL", "ERVK", "ERV1"):
        return "LTR"
    return "OTHER"


def is_erv(repeat_class: str) -> bool:
    """True if the locus is an endogenous-retrovirus (ERV) LTR element.

    ERVs are the LTR subset with retroviral gag/pol/env coding potential —
    the family the task flags as the priority antigen source. We require the
    locus to be in the LTR bucket AND carry an ERV/HERV-family token.
    """
    if repeat_class is None:
        return False
    s = str(repeat_class).strip().upper()
    if classify_family(s) != "LTR":
        return False
    return any(tok in s for tok in _ERV_TOKENS)


# ---------------------------------------------------------------------------
# 1) EXPRESSED / ACTIVE LOCUS SELECTION  (within-sample, multimap-aware)
# ---------------------------------------------------------------------------
def select_expressed_loci(counts: Mapping[str, float] | pd.Series,
                          min_reads: float = 10.0,
                          min_cpm: float = 1.0) -> list[str]:
    """Select transcriptionally ACTIVE TE loci for one sample.

    Telescope reports an EM-reassigned read count per TE locus. A locus is
    called active if it clears BOTH an absolute floor (``min_reads``, to reject
    single-read multimap noise) AND a WITHIN-sample CPM floor (``min_cpm``,
    library-size normalized so the call does not drift with sequencing depth).

    Within-sample CPM = 1e6 * count / (sum of all TE-locus counts in the
    sample). This is deliberately self-normalizing: it is the batch-robust way
    to call activity because it does not depend on absolute depth or on
    gene-level library size (the project's reproducibility constraint).

    Parameters
    ----------
    counts   : {locus_id: reassigned_read_count} for ONE sample (Telescope).
    min_reads: absolute read floor (default 10) — rejects multimap noise.
    min_cpm  : within-sample CPM floor (default 1.0).

    Returns the locus ids that pass, sorted by descending count.
    """
    s = pd.Series(dict(counts), dtype="float64") if not isinstance(counts, pd.Series) \
        else counts.astype("float64")
    s = s[s > 0]
    if s.empty:
        return []
    total = float(s.sum())
    cpm = 1e6 * s / total if total > 0 else s * 0.0
    keep = s[(s >= min_reads) & (cpm >= min_cpm)]
    return list(keep.sort_values(ascending=False).index)


# ---------------------------------------------------------------------------
# 2) ORF -> PEPTIDE DERIVATION  (frame translation of expressed locus sequence)
# ---------------------------------------------------------------------------
_CODON_TABLE = 1  # standard genetic code (Bio.Data.CodonTable id 1)


def _translate_frames(seq: str, strand: str = ".") -> list[str]:
    """Translate a nucleotide sequence in the relevant reading frames.

    strand '+' -> 3 forward frames; '-' -> 3 reverse-complement frames;
    '.' or unknown -> all 6 frames (we don't assume the element's orientation).
    Returns the list of raw protein strings (may contain '*' stop symbols).
    """
    from Bio.Seq import Seq

    s = "".join(c for c in str(seq).strip().upper() if c in "ACGTN")
    if len(s) < 3 * MIN_PEPTIDE_LEN:
        return []
    proteins: list[str] = []
    strands: list[str]
    if strand == "+":
        strands = ["+"]
    elif strand == "-":
        strands = ["-"]
    else:
        strands = ["+", "-"]
    for st in strands:
        nt = Seq(s) if st == "+" else Seq(s).reverse_complement()
        for frame in range(3):
            sub = nt[frame:]
            sub = sub[: len(sub) - (len(sub) % 3)]
            if len(sub) < 3 * MIN_PEPTIDE_LEN:
                continue
            # translate with '*' at stops; do not stop early
            proteins.append(str(sub.translate(table=_CODON_TABLE, to_stop=False)))
    return proteins


def _orf_segments(protein: str, min_orf_aa: int) -> list[str]:
    """Split a translated frame on stop codons into candidate ORF segments.

    Each maximal run of coding residues between stops is a candidate ORF. We
    keep segments at least ``min_orf_aa`` long. We do NOT require an internal
    Met start: TE/ERV ORFs frequently initiate at non-canonical or upstream
    starts, and MHC-I peptides are generated by proteasomal cleavage of the
    translated product regardless of the initiator — requiring 'M' would drop
    real antigens. (A Met-restricted variant is available via require_met.)
    """
    out = []
    for seg in protein.split("*"):
        if len(seg) >= min_orf_aa:
            out.append(seg)
    return out


def peptides_from_sequence(seq: str,
                           strand: str = ".",
                           kmin: int = MIN_PEPTIDE_LEN,
                           kmax: int = MAX_PEPTIDE_LEN,
                           min_orf_aa: int | None = None,
                           require_met: bool = False) -> set[str]:
    """All candidate 8-11mer peptides from the ORFs of one locus sequence.

    Steps: translate the locus sequence in the relevant frame(s) ->
    split each frame on stops into ORF segments (>= min_orf_aa) ->
    (optionally trim each ORF to its first Met) -> tile every kmin..kmax-mer.

    ``min_orf_aa`` defaults to ``kmin`` (an ORF must be at least one peptide
    long). Peptide hygiene (standard-AA, dedup, length) is enforced again by
    the shared engine, so returning a superset here is safe.
    """
    if min_orf_aa is None:
        min_orf_aa = kmin
    peps: set[str] = set()
    for protein in _translate_frames(seq, strand=strand):
        for orf in _orf_segments(protein, min_orf_aa):
            if require_met:
                i = orf.find("M")
                if i < 0:
                    continue
                orf = orf[i:]
                if len(orf) < min_orf_aa:
                    continue
            L = len(orf)
            for k in range(kmin, kmax + 1):
                if L < k:
                    break
                for i in range(0, L - k + 1):
                    kmer = orf[i:i + k]
                    if "X" not in kmer:          # engine drops it anyway
                        peps.add(kmer)
    return peps


def peptides_by_locus(locus_seqs: Mapping[str, str],
                      strand_map: Mapping[str, str] | None = None,
                      **kwargs) -> dict[str, set[str]]:
    """Map each expressed locus id -> its candidate peptide set.

    locus_seqs : {locus_id: genomic_sequence}. strand_map (optional):
    {locus_id: '+'/'-'} from the annotation; loci absent from it are
    translated in all 6 frames.
    """
    strand_map = strand_map or {}
    out: dict[str, set[str]] = {}
    for lid, seq in locus_seqs.items():
        out[lid] = peptides_from_sequence(
            seq, strand=strand_map.get(lid, "."), **kwargs)
    return out


# ---------------------------------------------------------------------------
# 3) GENOMIC SEQUENCE EXTRACTION for expressed loci (optional convenience)
# ---------------------------------------------------------------------------
def extract_locus_sequences(annotation: pd.DataFrame,
                            genome_fasta: str | Path,
                            locus_ids: Sequence[str] | None = None,
                            id_col: str = "locus_id",
                            chrom_col: str = "chrom",
                            start_col: str = "start",
                            end_col: str = "end") -> dict[str, str]:
    """Pull genomic sequence for each locus from a genome FASTA (pysam).

    ``annotation`` must have id/chrom/start/end columns (0-based half-open
    start, as in BED/RepeatMasker-to-BED). Only ``locus_ids`` are extracted if
    given. Returns {locus_id: sequence}. Requires the genome FASTA + .fai; on
    the pilot host this is the same reference the pipeline aligned against.
    """
    import pysam

    ann = annotation
    if locus_ids is not None:
        ann = ann[ann[id_col].isin(set(locus_ids))]
    seqs: dict[str, str] = {}
    fa = pysam.FastaFile(str(genome_fasta))
    try:
        for _, r in ann.iterrows():
            chrom = str(r[chrom_col])
            start = int(r[start_col]); end = int(r[end_col])
            try:
                seqs[str(r[id_col])] = fa.fetch(chrom, start, end)
            except (KeyError, ValueError):
                # contig not in this FASTA build — skip, don't fabricate
                continue
    finally:
        fa.close()
    return seqs


# ---------------------------------------------------------------------------
# 4) THE BURDEN — per sample, via the SHARED engine
# ---------------------------------------------------------------------------
_ROW_COLS = [
    "run_accession", "cohort",
    "te_antigen_burden", "te_antigen_burden_strong",
    "te_antigen_burden_LINE", "te_antigen_burden_SINE",
    "te_antigen_burden_LTR", "te_antigen_burden_ERV",
    "te_antigen_n_expressed_loci", "te_antigen_n_binder_loci",
    "te_antigen_top_locus",
]


def compute_te_antigen_burden(
        counts: Mapping[str, float] | pd.Series,
        locus_seqs: Mapping[str, str],
        hla_alleles: Sequence[str],
        annotation: pd.DataFrame | None = None,
        family_map: Mapping[str, str] | None = None,
        strand_map: Mapping[str, str] | None = None,
        min_reads: float = 10.0,
        min_cpm: float = 1.0,
        rank_threshold: float = WEAK_BINDER_RANK,
        id_col: str = "locus_id",
        class_col: str = "repeat_class",
        strand_col: str = "strand") -> dict:
    """Per-sample TE/ERV antigen burden through the shared MHC-I engine.

    Parameters
    ----------
    counts      : {locus_id: Telescope reassigned count} for one sample.
    locus_seqs  : {locus_id: genomic sequence}. Only expressed loci need a
                  sequence; extra sequences are ignored.
    hla_alleles : the sample's HLA-I alleles (from antigen_core.hla_typing).
    annotation  : optional DataFrame with `id_col`, `class_col`, `strand_col`
                  used to derive family + strand maps if those are not passed.
    family_map  : optional {locus_id: repeat_class} (overrides annotation).
    strand_map  : optional {locus_id: '+'/'-'}.

    Returns a dict with the NAMED feature columns (see _ROW_COLS) plus:
        'locus_contributions' : DataFrame [locus_id, family, is_erv, count,
                                 n_peptides, n_binders] for expressed loci,
        'scored'              : the per-peptide DataFrame from the engine.

    Method
    ------
    active loci (within-sample CPM filter) -> ORF peptides per locus ->
    POOL unique peptides -> score ONCE against the sample's HLA genotype ->
    a peptide is a binder if affinity percentile <= rank_threshold. The
    headline burden counts UNIQUE binder peptides across all active loci
    (== antigen_core.count_binders on the pooled set, asserted in tests).
    Family/locus burdens attribute each binder peptide to the loci that
    produced it (a peptide from multiple loci is counted once per family it
    appears in, and once in the headline).
    """
    # --- resolve family + strand maps -------------------------------------
    if family_map is None:
        family_map = {}
        if annotation is not None and class_col in annotation.columns:
            family_map = dict(zip(annotation[id_col], annotation[class_col]))
    if strand_map is None:
        strand_map = {}
        if annotation is not None and strand_col in annotation.columns:
            strand_map = dict(zip(annotation[id_col], annotation[strand_col]))

    # --- 1) active loci ----------------------------------------------------
    active = select_expressed_loci(counts, min_reads=min_reads, min_cpm=min_cpm)
    counts_s = pd.Series(dict(counts), dtype="float64")

    # --- 2) peptides per active locus (only those with sequence) ----------
    active_with_seq = [l for l in active if l in locus_seqs and locus_seqs[l]]
    pep_by_locus = peptides_by_locus(
        {l: locus_seqs[l] for l in active_with_seq},
        strand_map=strand_map,
    )

    # pooled unique candidate peptides
    all_peps: set[str] = set()
    for ps in pep_by_locus.values():
        all_peps |= ps

    # --- 3) score once via the SHARED engine ------------------------------
    scored = best_per_peptide(sorted(all_peps), list(hla_alleles))
    if scored.empty:
        binder_peps: set[str] = set()
    else:
        binder_peps = set(
            scored.loc[scored["affinity_percentile"] <= rank_threshold, "peptide"]
        )
    strong_peps = set() if scored.empty else set(
        scored.loc[scored["affinity_percentile"] <= STRONG_BINDER_RANK, "peptide"]
    )

    # --- 4) attribute binders to loci + families --------------------------
    fam_binders: dict[str, set[str]] = {f: set() for f in FAMILY_BUCKETS}
    erv_binders: set[str] = set()
    rows = []
    for lid in active_with_seq:
        peps = pep_by_locus.get(lid, set())
        b = peps & binder_peps
        fam = classify_family(family_map.get(lid, ""))
        erv = is_erv(family_map.get(lid, ""))
        fam_binders[fam] |= b
        if erv:
            erv_binders |= b
        rows.append({
            "locus_id": lid,
            "family": fam,
            "is_erv": erv,
            "count": float(counts_s.get(lid, 0.0)),
            "n_peptides": len(peps),
            "n_binders": len(b),
        })
    contrib = pd.DataFrame(rows, columns=[
        "locus_id", "family", "is_erv", "count", "n_peptides", "n_binders"])
    if not contrib.empty:
        contrib = contrib.sort_values(
            ["n_binders", "count"], ascending=False).reset_index(drop=True)

    n_binder_loci = int((contrib["n_binders"] > 0).sum()) if not contrib.empty else 0
    top_locus = (contrib.iloc[0]["locus_id"]
                 if (not contrib.empty and contrib.iloc[0]["n_binders"] > 0)
                 else "")

    result = {
        "run_accession": None,
        "cohort": None,
        "te_antigen_burden": len(binder_peps),
        "te_antigen_burden_strong": len(strong_peps),
        "te_antigen_burden_LINE": len(fam_binders["LINE"]),
        "te_antigen_burden_SINE": len(fam_binders["SINE"]),
        "te_antigen_burden_LTR": len(fam_binders["LTR"]),
        "te_antigen_burden_ERV": len(erv_binders),
        "te_antigen_n_expressed_loci": len(active),
        "te_antigen_n_binder_loci": n_binder_loci,
        "te_antigen_top_locus": top_locus,
        "locus_contributions": contrib,
        "scored": scored,
    }
    return result


def te_antigen_row(run_accession: str,
                   cohort: str,
                   counts: Mapping[str, float] | pd.Series,
                   locus_seqs: Mapping[str, str],
                   hla_alleles: Sequence[str],
                   **kwargs) -> dict:
    """One tidy feature-matrix row (keyed run_accession+cohort).

    Thin wrapper over ``compute_te_antigen_burden`` that stamps the sample key
    and returns only the scalar named columns (_ROW_COLS) — the shape written
    into results/features/. The locus_contributions / scored tables are dropped
    (kept only in the full compute result for provenance/QC).
    """
    r = compute_te_antigen_burden(counts, locus_seqs, hla_alleles, **kwargs)
    r["run_accession"] = run_accession
    r["cohort"] = cohort
    return {c: r[c] for c in _ROW_COLS}


def build_te_antigen_table(rows: Iterable[dict]) -> pd.DataFrame:
    """Assemble per-sample rows into the contract-format tidy matrix."""
    df = pd.DataFrame(list(rows))
    for c in _ROW_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    return df[_ROW_COLS].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Telescope report parser (provenance helper — the pipeline session emits
# te_locus.parquet, but this lets the module ingest a raw Telescope report
# directly for validation / one-off runs).
# ---------------------------------------------------------------------------
def parse_telescope_report(path: str | Path,
                           count_col: str = "final_count") -> pd.Series:
    """Parse a Telescope ``*-telescope_report.tsv`` -> {locus_id: count}.

    Telescope writes a '## RunInfo' comment line, then a header row with
    columns: transcript, transcript_length, final_count, final_conf,
    final_prop, init_aligned, unique_count, init_best, ... We keep
    ``transcript`` as the locus id and ``final_count`` (EM-reassigned reads).
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    idc = "transcript" if "transcript" in df.columns else df.columns[0]
    if count_col not in df.columns:
        # fall back to first integer count-like column
        count_col = next((c for c in df.columns if "count" in c.lower()), df.columns[2])
    s = pd.Series(df[count_col].values, index=df[idc].astype(str).values,
                  dtype="float64")
    return s[s.index != "__no_feature"]


# ---------------------------------------------------------------------------
# BATCH ROBUSTNESS NOTE (project hard constraint — built in from the start)
# ---------------------------------------------------------------------------
def build_batch_robustness_note() -> str:
    """Return the documented batch/platform-reproducibility note for this feature."""
    return (
        "TE/ERV quantification is MULTIMAP-SENSITIVE: TE reads map to many "
        "near-identical genomic copies, so the per-locus count depends entirely "
        "on the aligner's multimapping policy. REQUIREMENT: every sample's "
        "Telescope (or TEtranscripts) counts MUST come from the SAME STAR "
        "multimap settings — specifically the same "
        "--outFilterMultimapNmax / --winAnchorMultimapNmax / "
        "--outSAMmultNmax and Telescope --theta_prior across the whole cohort; "
        "mixing settings makes locus counts non-comparable. "
        "MITIGATIONS BUILT INTO THIS MODULE: (1) the activity call is a "
        "WITHIN-sample CPM normalization (1e6*count/sum-of-TE-counts), so it "
        "does not drift with sequencing depth or gene-level library size; "
        "(2) an absolute min-read floor rejects single-read multimap noise; "
        "(3) the burden itself is an MHCflurry percentile-RANK binder count "
        "(allele-calibrated against a fixed random-peptide background), which "
        "is composition- and platform-invariant. "
        "READ-LENGTH CAVEAT: TE quantification is read-length sensitive; "
        "validate the feature per platform and z-score within cohort before "
        "pooling — never pool raw counts across cohorts."
    )


__all__ = [
    "classify_family", "is_erv", "FAMILY_BUCKETS",
    "select_expressed_loci",
    "peptides_from_sequence", "peptides_by_locus", "extract_locus_sequences",
    "compute_te_antigen_burden", "te_antigen_row", "build_te_antigen_table",
    "parse_telescope_report", "build_batch_robustness_note",
]
