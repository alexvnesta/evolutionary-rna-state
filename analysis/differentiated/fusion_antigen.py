"""
analysis/differentiated/fusion_antigen.py

DIFFERENTIATED BUCKET — Fusion-transcript-derived neoantigen burden.

Gene fusions create chimeric transcripts whose junction-spanning reading frame
encodes a peptide sequence that exists in NO wild-type protein. The amino acids
straddling the fusion breakpoint are therefore tumour-specific neo-epitopes,
invisible to any SNV/indel-based (WES) neoantigen predictor. This module turns a
sample's fusion calls into a single, named, interpretable burden feature by
routing the junction peptides through the SAME shared MHCflurry engine every
other antigen module uses.

Pipeline (per sample)
---------------------
    Arriba fusions.tsv  (or STAR-Fusion coding-effect TSV)
      -> keep IN-FRAME fusions (reading_frame == 'in-frame' / PROT_FUSION_TYPE
         == 'INFRAME'), high/medium confidence
      -> take the translated fusion-junction ORF (Arriba `peptide_sequence`,
         breakpoint marked '|'; or translate `fusion_transcript`)
      -> enumerate 8-11mer peptides that CROSS the breakpoint (>=1 residue on
         each side of the junction) — these are the novel epitopes
      -> antigen_core.mhc_binding.count_binders(peptides, hla_alleles)
      -> fusion_neoantigen_burden  (int, per sample)

Named interpretable features emitted (v2 contract)
--------------------------------------------------
    fusion_neoantigen_burden        int   unique MHC-I binders across all
                                          in-frame fusion junctions (rank<=2.0)
    fusion_neoantigen_burden_strong int   strong-only (rank<=0.5), for QC
    n_fusions                       int   total fusion calls (passing confidence)
    n_inframe_fusions               int   subset that are protein-coding in-frame

Why the burden is comparable to the sibling burdens: it is a
``count_binders`` call on the ONE shared engine (MHCflurry percentile rank,
allele-calibrated), identical definition to splice/TE/IR/SNV burdens.

Batch-robustness (fusion calling is caller- and depth-sensitive)
----------------------------------------------------------------
See ``BATCH_ROBUSTNESS_NOTE`` at the bottom of this file. In short: FIX the
caller + version + reference across every sample in a comparison, filter to
high/medium confidence, and never mix Arriba-derived and STAR-Fusion-derived
burdens in the same axis. The MHC engine (percentile rank) is batch-invariant;
the *fusion detection* upstream is not, so the caller must be held constant.

References
----------
Uhrig et al. 2021, Genome Research, "Accurate and efficient detection of gene
fusions from RNA sequencing data" (Arriba; doi:10.1101/gr.257246.119).
Haas et al. 2019, Genome Biology, "Accuracy assessment of fusion transcript
detection via read-mapping and de novo fusion transcript assembly-based
methods" (STAR-Fusion; doi:10.1186/s13059-019-1842-9).
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Import the SHARED antigen engine. Works both as a package
# (analysis.differentiated.fusion_antigen) and when the file's directory is on
# sys.path (Nextflow work dir / bin script), mirroring bin/merge_hla_table.py.
# ---------------------------------------------------------------------------
try:
    from analysis.antigen_core.mhc_binding import (
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
except ImportError:  # pragma: no cover - path shim for Nextflow / standalone
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from antigen_core.mhc_binding import (  # type: ignore
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
K_MIN = 8          # class-I peptide window (matches the shared engine)
K_MAX = 11
DEFAULT_CONFIDENCE = ("high", "medium")   # Arriba confidence tiers to keep
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# reading-frame labels that mean "protein-coding, in-frame across the junction"
_INFRAME_ARRIBA = {"in-frame"}
_INFRAME_STARFUSION = {"INFRAME"}

FEATURE_COLS = [
    "run_accession", "cohort",
    "n_fusions", "n_inframe_fusions",
    "fusion_neoantigen_burden", "fusion_neoantigen_burden_strong",
    "caller", "caller_version",
]


# ---------------------------------------------------------------------------
# Unified per-fusion record
# ---------------------------------------------------------------------------
@dataclass
class FusionCall:
    """One caller-agnostic fusion call.

    ``peptide`` is the translated junction-spanning protein with the breakpoint
    marked by '|'. ``transcript`` is the fused nucleotide sequence (also '|'-
    marked) — used to translate the ORF when a peptide was not pre-computed.
    """
    gene1: str
    gene2: str
    breakpoint1: str = ""
    breakpoint2: str = ""
    reading_frame: str = ""          # caller's raw label
    confidence: str = ""
    peptide: str = ""                # '|'-marked translated junction peptide
    transcript: str = ""             # '|'-marked fused nucleotide sequence
    support: int = 0                 # supporting reads (split + discordant)
    caller: str = ""

    @property
    def name(self) -> str:
        return f"{self.gene1}--{self.gene2}"

    def is_inframe(self) -> bool:
        rf = str(self.reading_frame).strip()
        return rf in _INFRAME_ARRIBA or rf in _INFRAME_STARFUSION


# ---------------------------------------------------------------------------
# Peptide extraction around the breakpoint
# ---------------------------------------------------------------------------
def _window_peptides(protein: str, junction_idx: int,
                     k_min: int = K_MIN, k_max: int = K_MAX) -> list[str]:
    """All k-mers (k in [k_min,k_max]) that CROSS ``junction_idx``.

    ``junction_idx`` is the number of residues to the LEFT of the breakpoint,
    i.e. residues protein[0:junction_idx] come from gene1 and
    protein[junction_idx:] from gene2. A window [start,end) crosses the junction
    iff it contains >=1 residue from each side: start < junction_idx < end.

    Non-standard characters (stop '*', ambiguous '?', frameshift markers) are
    left in place here and dropped downstream by the engine's ``clean_peptides``
    (which enforces the 8-11 length and standard-AA alphabet), so a window that
    straddles a stop codon is simply not counted — the correct behaviour.
    """
    prot = protein
    L = len(prot)
    out: list[str] = []
    seen: set[str] = set()
    for k in range(k_min, k_max + 1):
        for start in range(0, L - k + 1):
            end = start + k
            if start < junction_idx < end:
                pep = prot[start:end]
                if pep not in seen:
                    seen.add(pep)
                    out.append(pep)
    return out


def peptides_from_marked_peptide(marked: str,
                                 k_min: int = K_MIN,
                                 k_max: int = K_MAX) -> list[str]:
    """Junction-crossing k-mers from a '|'-marked translated peptide.

    Handles the Arriba ``peptide_sequence`` grammar:
      '|'   fusion breakpoint (required; the residue AFTER '|' is the first
            gene2 residue)
      '...' ORF extends beyond the shown window -> stripped
      '*'   stop codon -> the peptide is truncated here (windows past the stop
            are not generated)
      lower-case (frameshift-region residues) -> upper-cased; these ARE real
            translated residues and are the most tumour-specific epitopes
      '?'   ambiguous residue at a split codon -> left as-is, so any window
            containing it is dropped by the engine (conservative)

    Returns [] when there is no breakpoint marker or no coding sequence.
    """
    if not marked or str(marked).strip() in (".", ""):
        return []
    s = str(marked).strip()
    if "|" not in s:
        return []
    # first breakpoint marker defines the junction
    left, right = s.split("|", 1)
    # strip truncation ellipses and stray dots/whitespace, keep letter markers
    clean = lambda x: x.replace("...", "").replace(".", "").strip()
    left, right = clean(left), clean(right)
    full = (left + right).upper()
    junction = len(left)
    if junction == 0 or junction >= len(full):
        return []
    # truncate at first stop codon; if the stop is at/left of the junction there
    # is no viable crossing peptide
    if "*" in full:
        stop = full.index("*")
        if stop <= junction:
            return []
        full = full[:stop]
    return _window_peptides(full, junction, k_min, k_max)


def translate_fusion_transcript(transcript: str,
                                frame: int | None = None) -> tuple[str, int]:
    """Translate a '|'-marked fused nucleotide transcript to (protein, junction).

    Used when a caller supplies only a nucleotide fusion transcript (no pre-
    translated peptide). The '|' marks the breakpoint in nucleotide space. We:
      1. record the nucleotide junction index (bases left of '|'),
      2. pick the reading frame: if ``frame`` is given use it; otherwise choose
         the frame (0/1/2) that yields the longest ORF with no premature stop
         before the junction (falls back to frame 0),
      3. translate with Biopython,
      4. map the nucleotide junction to an amino-acid index (nt_junction//3 in
         the chosen frame).

    Returns (protein, aa_junction_index). Lower-case / marker characters in the
    input transcript are stripped; only A/C/G/T/N are translated.
    """
    if not transcript or str(transcript).strip() in (".", ""):
        return "", 0
    s = str(transcript).strip().upper()
    s = s.replace("...", "")
    nt_junction = s.index("|") if "|" in s else len(s)
    # keep only nucleotides; track how many were removed before the junction so
    # the junction index stays correct
    left_raw = s[:nt_junction]
    left_nt = re.sub(r"[^ACGTN]", "", left_raw)
    all_nt = re.sub(r"[^ACGTN]", "", s)
    nt_junction = len(left_nt)

    def _translate(f: int) -> str:
        seq = all_nt[f:]
        seq = seq[: len(seq) - (len(seq) % 3)]
        if not seq:
            return ""
        return str(Seq(seq).translate())

    if frame is None:
        best_f, best_len = 0, -1
        for f in (0, 1, 2):
            prot = _translate(f)
            aa_j = (nt_junction - f) // 3
            if aa_j < 0:
                continue
            pre = prot[:aa_j]
            score = len(prot) if "*" not in pre else aa_j  # penalise early stop
            if score > best_len:
                best_f, best_len = f, score
        frame = best_f
    prot = _translate(frame)
    aa_junction = max((nt_junction - frame) // 3, 0)
    return prot, aa_junction


def fusion_peptides(call: FusionCall,
                    k_min: int = K_MIN, k_max: int = K_MAX) -> list[str]:
    """Junction-crossing peptides for one IN-FRAME fusion call.

    Prefers the caller's pre-translated ``peptide`` (Arriba peptide_sequence /
    STAR-Fusion FUSION_TRANSLATION with a '|' junction). Falls back to
    translating the nucleotide ``transcript`` ORF. Returns [] for out-of-frame
    calls or when no coding junction sequence is available.
    """
    if not call.is_inframe():
        return []
    if call.peptide and "|" in str(call.peptide):
        peps = peptides_from_marked_peptide(call.peptide, k_min, k_max)
        if peps:
            return peps
    if call.transcript:
        prot, junc = translate_fusion_transcript(call.transcript)
        if prot and 0 < junc < len(prot):
            # stop translation at the first stop codon at/after the junction
            if "*" in prot[junc:]:
                stop = junc + prot[junc:].index("*")
                prot = prot[:stop]
            if junc < len(prot):
                return _window_peptides(prot, junc, k_min, k_max)
    return []


# ---------------------------------------------------------------------------
# Caller parsers -> list[FusionCall]
# ---------------------------------------------------------------------------
def _norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def read_arriba(path: str | Path,
                confidence_levels: Sequence[str] = DEFAULT_CONFIDENCE
                ) -> list[FusionCall]:
    """Parse an Arriba ``fusions.tsv`` into FusionCall records.

    Arriba columns used: gene1, gene2, breakpoint1, breakpoint2, confidence,
    reading_frame, peptide_sequence, fusion_transcript, split_reads1,
    split_reads2, discordant_mates. Rows are kept only if confidence is in
    ``confidence_levels`` (default high+medium) — the batch-robust filter.
    """
    df = _norm_cols(pd.read_csv(path, sep="\t", dtype=str).fillna(""))
    keep = {c.lower(): c for c in df.columns}
    def col(name, default=""):
        return df[keep[name]] if name in keep else pd.Series([default] * len(df))
    calls: list[FusionCall] = []
    conf_set = {c.lower() for c in confidence_levels}
    for _, r in df.iterrows():
        conf = str(r.get(keep.get("confidence", ""), "")).strip().lower()
        if conf_set and conf not in conf_set:
            continue
        def g(n):
            return str(r[keep[n]]).strip() if n in keep else ""
        try:
            support = (int(float(g("split_reads1") or 0))
                       + int(float(g("split_reads2") or 0))
                       + int(float(g("discordant_mates") or 0)))
        except ValueError:
            support = 0
        calls.append(FusionCall(
            gene1=g("gene1"), gene2=g("gene2"),
            breakpoint1=g("breakpoint1"), breakpoint2=g("breakpoint2"),
            reading_frame=g("reading_frame"),
            confidence=conf,
            peptide=g("peptide_sequence"),
            transcript=g("fusion_transcript"),
            support=support,
            caller="arriba",
        ))
    return calls


def read_starfusion(path: str | Path) -> list[FusionCall]:
    """Parse a STAR-Fusion coding-effect TSV into FusionCall records.

    Requires STAR-Fusion run with ``--examine_coding_effect``. Columns used:
    #FusionName, LeftBreakpoint, RightBreakpoint, PROT_FUSION_TYPE
    (INFRAME/FRAMESHIFT), FUSION_TRANSLATION, CDS_LEFT_RANGE (to locate the
    junction), JunctionReadCount, SpanningFragCount.

    STAR-Fusion's FUSION_TRANSLATION is not '|'-marked, so the junction is
    reconstructed from the left CDS length (aa_junction = left_cds_len // 3) and
    a '|' is inserted so downstream peptide extraction is caller-agnostic.
    """
    df = _norm_cols(pd.read_csv(path, sep="\t", dtype=str).fillna(""))
    keep = {c.lower(): c for c in df.columns}
    calls: list[FusionCall] = []
    for _, r in df.iterrows():
        def g(n):
            return str(r[keep[n]]).strip() if n in keep else ""
        name = g("fusionname")
        g1, _, g2 = name.partition("--")
        prot = g("fusion_translation")
        # locate junction from CDS_LEFT_RANGE ("start-end", 1-based on the CDS)
        marked = prot
        cds_left = g("cds_left_range")
        m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", cds_left)
        if prot and prot != "." and m:
            left_len = int(m.group(2)) - int(m.group(1)) + 1
            aa_j = left_len // 3
            if 0 < aa_j < len(prot):
                marked = prot[:aa_j] + "|" + prot[aa_j:]
        try:
            support = (int(float(g("junctionreadcount") or 0))
                       + int(float(g("spanningfragcount") or 0)))
        except ValueError:
            support = 0
        calls.append(FusionCall(
            gene1=g1, gene2=g2,
            breakpoint1=g("leftbreakpoint"), breakpoint2=g("rightbreakpoint"),
            reading_frame=g("prot_fusion_type"),
            confidence="", peptide=marked, transcript="",
            support=support, caller="starfusion",
        ))
    return calls


# ---------------------------------------------------------------------------
# Per-sample feature computation
# ---------------------------------------------------------------------------
def compute_fusion_features(calls: Sequence[FusionCall],
                            hla_alleles: Sequence[str],
                            k_min: int = K_MIN,
                            k_max: int = K_MAX) -> dict:
    """Core: FusionCall list + HLA alleles -> named feature dict.

    Pools the junction peptides from ALL in-frame fusions and scores them once
    through the shared engine, so ``fusion_neoantigen_burden`` counts UNIQUE
    binding peptides across the sample (a peptide arising from two fusions is
    one antigen — consistent with the sibling burdens).
    """
    n_fusions = len(calls)
    inframe = [c for c in calls if c.is_inframe()]
    n_inframe = len(inframe)

    peptides: list[str] = []
    for c in inframe:
        peptides.extend(fusion_peptides(c, k_min, k_max))
    # dedupe while preserving order (engine also dedupes/cleans)
    peptides = list(dict.fromkeys(peptides))

    if peptides and len(hla_alleles) > 0:
        burden = count_binders(peptides, hla_alleles, rank_threshold=WEAK_BINDER_RANK)
        burden_strong = count_binders(peptides, hla_alleles, rank_threshold=STRONG_BINDER_RANK)
    else:
        burden = 0
        burden_strong = 0

    return {
        "n_fusions": int(n_fusions),
        "n_inframe_fusions": int(n_inframe),
        "fusion_neoantigen_burden": int(burden),
        "fusion_neoantigen_burden_strong": int(burden_strong),
        "n_junction_peptides": int(len(peptides)),  # QC (not a contract column)
    }


def fusion_features_for_sample(run_accession: str,
                               cohort: str,
                               hla_alleles: Sequence[str],
                               arriba_tsv: str | Path | None = None,
                               starfusion_tsv: str | Path | None = None,
                               confidence_levels: Sequence[str] = DEFAULT_CONFIDENCE,
                               caller_version: str = "") -> dict:
    """End-to-end for one sample -> a v2-contract feature row.

    Exactly ONE caller input must be given (the fixed-caller requirement — see
    the batch-robustness note). Returns a dict with FEATURE_COLS.
    """
    if bool(arriba_tsv) == bool(starfusion_tsv):
        raise ValueError(
            "Provide exactly one caller input (arriba_tsv XOR starfusion_tsv). "
            "Mixing callers across samples breaks batch comparability — see "
            "BATCH_ROBUSTNESS_NOTE."
        )
    if arriba_tsv:
        calls = read_arriba(arriba_tsv, confidence_levels=confidence_levels)
        caller = "arriba"
    else:
        calls = read_starfusion(starfusion_tsv)
        caller = "starfusion"

    feats = compute_fusion_features(calls, hla_alleles)
    row = {
        "run_accession": run_accession,
        "cohort": cohort,
        "n_fusions": feats["n_fusions"],
        "n_inframe_fusions": feats["n_inframe_fusions"],
        "fusion_neoantigen_burden": feats["fusion_neoantigen_burden"],
        "fusion_neoantigen_burden_strong": feats["fusion_neoantigen_burden_strong"],
        "caller": caller,
        "caller_version": caller_version,
    }
    return row


def build_fusion_feature_table(rows: Iterable[dict]) -> pd.DataFrame:
    """Assemble per-sample rows into the tidy v2-contract table."""
    df = pd.DataFrame(list(rows))
    for c in FEATURE_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    return df[FEATURE_COLS].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Batch-robustness note (fusion calling is caller/depth sensitive)
# ---------------------------------------------------------------------------
BATCH_ROBUSTNESS_NOTE = """\
FUSION_NEOANTIGEN_BURDEN — batch/platform robustness
====================================================
The MHC-binding half of this feature is batch-invariant (MHCflurry percentile
rank is calibrated per allele against a fixed random-peptide background, so it
does not drift with input composition). The FUSION-DETECTION half is NOT: fusion
call sets depend strongly on the caller, its version, the reference/annotation,
and — critically — sequencing depth and read length (more reads -> more junction
support -> more low-confidence calls). Uncontrolled, this injects a technical
axis that mimics biological signal.

FIXED-CALLER REQUIREMENT (enforced in code):
  * `fusion_features_for_sample` accepts exactly ONE caller input (Arriba XOR
    STAR-Fusion) and records `caller` + `caller_version` on every row. Never
    mix Arriba-derived and STAR-Fusion-derived burdens on the same comparison
    axis; the two callers have different sensitivity/precision profiles and
    different in-frame annotation logic.
  * Hold the caller VERSION and the reference/annotation build constant across
    every sample in a comparison (record both).

CONFIDENCE FILTER: default to Arriba high+medium confidence (drops the depth-
sensitive low-confidence tail). This is the main lever against depth-driven
false positives. STAR-Fusion analogue: require FFPM / junction+spanning read
support thresholds and `--examine_coding_effect`.

DEPTH: where library sizes differ substantially across a cohort, either
down-sample to a common depth before calling, or include per-sample fusion read
depth as a covariate. Report `n_fusions` and `n_inframe_fusions` alongside the
burden so a depth confound is visible (a cohort with systematically higher
n_fusions likely differs in depth, not biology).

REPORTING: z-score / rank within cohort before pooling across cohorts, exactly
as for the sibling burdens; never pool raw counts across platforms.
"""


if __name__ == "__main__":  # pragma: no cover
    print(BATCH_ROBUSTNESS_NOTE)
