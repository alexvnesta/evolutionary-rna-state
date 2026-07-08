"""
analysis/antigen_core/mhc_binding.py

SHARED ANTIGEN CORE — MHC class I binding/presentation engine.

This is THE peptide-scoring engine every downstream antigen module imports
(splicing, TE/ERV, intron-retention, fusion, SNV/indel). It wraps MHCflurry
2.x (Class1 presentation predictor: an affinity predictor + an antigen-
processing predictor combined into a presentation score) behind a small,
stable, documented interface so the siblings never touch MHCflurry directly.

Design contract (keep stable — siblings depend on it):
    score_peptides(peptides, hla_alleles)   -> tidy per (peptide x allele) DataFrame
    count_binders(peptides, hla_alleles, ...) -> int  (per-sample binder count)
    best_per_peptide(...)                    -> best allele per peptide (collapse)

Why MHCflurry (not netMHCpan): MHCflurry is open-source (Apache-2.0),
pip-installable, CPU-only, and returns calibrated percentile ranks plus a
presentation score in one call — no license wall, reproducible across
machines. Rank semantics follow the field convention:
    strong binder : percentile rank <= 0.5
    weak binder   : percentile rank <= 2.0
(lower rank = stronger binding; rank is the batch/platform-robust unit because
it is calibrated against a fixed set of random peptides per allele, so it does
not drift with input composition the way raw nM affinity does).

The predictor models are downloaded once with:
    mhcflurry-downloads fetch models_class1_presentation

References
---------
O'Donnell et al. 2020, Cell Systems, "MHCflurry 2.0: Improved Pan-Allele
Prediction of MHC Class I-Presented Peptides by Incorporating Antigen
Processing" (doi:10.1016/j.cels.2020.06.010).
"""
from __future__ import annotations

import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

# ---------------------------------------------------------------------------
# Model location. MHCflurry looks up its models via $MHCFLURRY_DATA_DIR. The
# macOS default (~/Library/Application Support/mhcflurry) is not writable in
# the sandbox, so the models were fetched into the repo's reference/ dir. If
# the caller has not set MHCFLURRY_DATA_DIR, point it there automatically so
# siblings importing this engine "just work" without extra setup.
# ---------------------------------------------------------------------------
_REPO_MODEL_DIR = (Path(__file__).resolve().parents[2]
                   / "reference" / "mhcflurry_models")
if "MHCFLURRY_DATA_DIR" not in os.environ and _REPO_MODEL_DIR.exists():
    os.environ["MHCFLURRY_DATA_DIR"] = str(_REPO_MODEL_DIR)

# ---------------------------------------------------------------------------
# Constants — the field-standard rank thresholds. Exposed so callers read them
# from here rather than hard-coding, keeping every module's binder definition
# identical.
# ---------------------------------------------------------------------------
STRONG_BINDER_RANK = 0.5      # percentile rank <= 0.5  -> strong binder
WEAK_BINDER_RANK = 2.0        # percentile rank <= 2.0  -> weak (incl. strong)
MIN_PEPTIDE_LEN = 8
MAX_PEPTIDE_LEN = 11          # MHCflurry class I supports 8-15; 8-11 is the
                              # biologically dominant class-I window we score.
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


@lru_cache(maxsize=1)
def _load_predictor():
    """Load the MHCflurry Class1 presentation predictor once (cached).

    Import is done lazily inside the function so that merely importing this
    module (e.g. for its constants) does not pull in TensorFlow. Raises a
    clear, actionable error if the models were never fetched.
    """
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError as e:  # pragma: no cover - env guard
        raise ImportError(
            "mhcflurry is not installed. Install with:\n"
            "    pip install mhcflurry\n"
            "    mhcflurry-downloads fetch models_class1_presentation"
        ) from e
    try:
        return Class1PresentationPredictor.load()
    except Exception as e:  # pragma: no cover - env guard
        raise RuntimeError(
            "Could not load MHCflurry models. Fetch them once with:\n"
            "    mhcflurry-downloads fetch models_class1_presentation"
        ) from e


def supported_alleles() -> list[str]:
    """Return the list of HLA alleles the loaded predictor supports."""
    return list(_load_predictor().supported_alleles)


# ---------------------------------------------------------------------------
# Peptide / allele cleaning
# ---------------------------------------------------------------------------
def clean_peptides(peptides: Iterable[str]) -> list[str]:
    """Uppercase, strip, dedupe, and keep only valid 8-11mers.

    Drops peptides with non-standard amino acids (X, U, *, gaps) and anything
    outside the 8-11 length window. Order is preserved (first occurrence).
    """
    seen: set[str] = set()
    out: list[str] = []
    for p in peptides:
        if p is None:
            continue
        pep = str(p).strip().upper()
        if not (MIN_PEPTIDE_LEN <= len(pep) <= MAX_PEPTIDE_LEN):
            continue
        if not _AA_RE.match(pep):
            continue
        if pep in seen:
            continue
        seen.add(pep)
        out.append(pep)
    return out


def normalize_allele(allele: str) -> str:
    """Normalize an HLA allele string to MHCflurry's expected form.

    Accepts common variants ('A*02:01', 'HLA-A*02:01', 'A0201', 'A_02_01')
    and returns 'HLA-A*02:01'. Returns the input unchanged if it cannot be
    confidently parsed (the predictor will then reject it, surfacing the issue).
    """
    a = str(allele).strip().upper().replace("HLA-", "")
    # A0201 -> A*02:01 ; A_02_01 -> A*02:01
    m = re.match(r"^([ABC])[\*_]?(\d{2})[:_]?(\d{2})$", a)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    m = re.match(r"^([ABC])\*(\d{2}):(\d{2})$", a)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    return str(allele).strip()


def _valid_alleles(hla_alleles: Sequence[str]) -> list[str]:
    """Normalize, dedupe, and keep only alleles the predictor supports."""
    supported = set(supported_alleles())
    out: list[str] = []
    seen: set[str] = set()
    for a in hla_alleles:
        na = normalize_allele(a)
        if na in seen:
            continue
        seen.add(na)
        if na in supported:
            out.append(na)
    return out


# ---------------------------------------------------------------------------
# Core scoring
# ---------------------------------------------------------------------------
_SCHEMA = ["peptide", "allele", "length", "affinity_nM",
           "affinity_percentile", "processing_score", "presentation_score",
           "presentation_percentile", "is_strong", "is_weak"]


def score_peptides(peptides: Sequence[str],
                   hla_alleles: Sequence[str]) -> pd.DataFrame:
    """Score peptides against a sample's HLA-I alleles.

    Uses MHCflurry's Class1 presentation predictor, which — given the sample's
    full genotype — reports for each peptide the SINGLE best-presenting allele
    (the biologically relevant call: a peptide is presented if any of the
    sample's alleles present it). Result is therefore one row per peptide.

    Parameters
    ----------
    peptides : list of candidate peptide strings (8-11mers). Cleaned/deduped
        internally; non-standard-AA and out-of-range peptides are dropped.
    hla_alleles : the sample's HLA-I alleles (up to 6: A/B/C x2). Accepts
        loose formats; normalized internally. Unsupported alleles are dropped.

    Returns
    -------
    DataFrame, one row per peptide, columns:
        peptide, allele (best presenter), length,
        affinity_nM             -- predicted binding affinity (lower = stronger)
        affinity_percentile     -- calibrated %-rank of affinity (lower = stronger)
        processing_score        -- antigen-processing model score
        presentation_score      -- combined processing+affinity score in [0,1]
                                   (higher = more likely presented)
        presentation_percentile -- %-rank of the presentation score
        is_strong               -- affinity_percentile <= STRONG_BINDER_RANK
        is_weak                 -- affinity_percentile <= WEAK_BINDER_RANK

    Empty DataFrame (with the schema) if no valid peptides or no supported
    alleles remain after cleaning.
    """
    peps = clean_peptides(peptides)
    alleles = _valid_alleles(hla_alleles)
    if not peps or not alleles:
        return pd.DataFrame(columns=_SCHEMA)

    predictor = _load_predictor()
    # alleles={sample: [...]} => scored against the whole genotype; predict()
    # returns one row per peptide with the best-presenting allele.
    raw = predictor.predict(
        peptides=peps,
        alleles={"sample": alleles},
        include_affinity_percentile=True,
        verbose=0,
    )

    df = pd.DataFrame({
        "peptide": raw["peptide"].values,
        "allele": raw["best_allele"].values,
        "affinity_nM": raw["affinity"].values,
        "affinity_percentile": raw["affinity_percentile"].values,
        "processing_score": raw["processing_score"].values,
        "presentation_score": raw["presentation_score"].values,
        "presentation_percentile": raw["presentation_percentile"].values,
    })
    df["length"] = df["peptide"].str.len()
    df["is_strong"] = df["affinity_percentile"] <= STRONG_BINDER_RANK
    df["is_weak"] = df["affinity_percentile"] <= WEAK_BINDER_RANK
    return df[_SCHEMA].reset_index(drop=True)


def best_per_peptide(peptides: Sequence[str],
                     hla_alleles: Sequence[str]) -> pd.DataFrame:
    """One row per peptide with its best-presenting allele.

    MHCflurry's presentation call already collapses to the best allele per
    peptide, so this returns ``score_peptides`` unchanged. Kept as an explicit,
    named entry point so downstream burden code reads clearly and stays
    correct if the underlying per-allele API is ever swapped in.
    """
    return score_peptides(peptides, hla_alleles)


def count_binders(peptides: Sequence[str],
                  hla_alleles: Sequence[str],
                  rank_threshold: float = WEAK_BINDER_RANK) -> int:
    """Per-sample BINDER COUNT: number of UNIQUE peptides that bind.

    A peptide counts as a binder if its best allele's affinity percentile is
    <= ``rank_threshold``. This is the atomic "*_neoantigen_burden" quantity
    every antigen module reports per sample.

    Parameters
    ----------
    rank_threshold : percentile-rank cutoff. Defaults to WEAK_BINDER_RANK
        (2.0). Pass STRONG_BINDER_RANK (0.5) for the strong-only count.

    Notes
    -----
    Counts UNIQUE peptides (a peptide binding several of the sample's alleles
    is one antigen), which is the biologically meaningful burden unit.
    """
    df = best_per_peptide(peptides, hla_alleles)
    if df.empty:
        return 0
    return int((df["affinity_percentile"] <= rank_threshold).sum())


def binder_counts(peptides: Sequence[str],
                  hla_alleles: Sequence[str]) -> dict[str, int]:
    """Convenience: both strong and weak binder counts in one pass.

    Returns {"n_strong_binders": int, "n_weak_binders": int, "n_scored": int}.
    ``n_weak_binders`` is inclusive of strong (rank <= 2.0), matching the
    field convention where "binder" defaults to the weak threshold.
    """
    df = best_per_peptide(peptides, hla_alleles)
    if df.empty:
        return {"n_strong_binders": 0, "n_weak_binders": 0, "n_scored": 0}
    return {
        "n_strong_binders": int((df["affinity_percentile"] <= STRONG_BINDER_RANK).sum()),
        "n_weak_binders": int((df["affinity_percentile"] <= WEAK_BINDER_RANK).sum()),
        "n_scored": int(len(df)),
    }
