"""
analysis/baseline/gep_scores.py
================================

BASELINE BUCKET — proven expression-signature floor for ICB response.

Named, batch-robust gene-expression signature scores computed from a
gene x sample TPM matrix. These are the *proven* baseline features that the
project's RNA-state neoantigen-burden features must beat. Three named scores:

    gep_tcell_inflamed  -- Ayers 2017 18-gene T-cell-inflamed GEP
    ifng_score          -- Ayers 2017 6-gene IFN-gamma signature
    teff_tgfb_balance   -- Mariathasan 2018 T-effector minus TGF-beta balance

All are floats, one value per sample, keyed on (run_accession, cohort).

Design decisions (batch/platform reproducibility is a HARD project constraint)
------------------------------------------------------------------------------
1.  INPUT is gene TPM. We work on log2(TPM + 1). log space stabilises variance
    and matches how all three signatures were originally defined (log-scale
    expression).

2.  HOUSEKEEPING NORMALISATION. The published Ayers GEP normalises each sample
    against a panel of housekeeping genes (subtract the sample's mean
    housekeeping log-expression) BEFORE combining signature genes. This removes
    per-sample library-size / input-amount scale and is the first line of
    platform robustness. We implement it exactly (`housekeeping=` genes) and
    fall back gracefully to no HK-centering if none of the HK genes are present,
    with a logged warning.

3.  WITHIN-BATCH STANDARDISATION (`harmonize=`). The single most important
    cross-platform guard: after computing the raw signature value, z-score (or
    rank-normalise) it WITHIN each cohort/batch so that a score is only ever
    compared to other samples measured on the same platform. Exposed as a switch
    because for a single-cohort analysis you may want the raw score. Default is
    within-cohort z-score.

4.  WEIGHTING. Ayers' T-cell-inflamed GEP is a *weighted* mean. The exact
    NanoString assay coefficients are published only in a proprietary
    clinical-assay supplement we could not retrieve from a primary open source,
    so the weight vector `AYERS_GEP_WEIGHTS_PUBLISHED` below is reconstructed
    from widely-reproduced secondary sources and is flagged UNVERIFIED. The
    DEFAULT weighting is therefore the equal-weight mean (`weighting="equal"`),
    which is fully reproducible and is the standard fallback used across the
    literature when the coefficients are unavailable. Pass
    `weighting="published"` to use the reconstructed coefficients. The gene LIST
    itself is verified (see provenance below).

Provenance
----------
* 18-gene T-cell-inflamed GEP & 6-gene IFN-gamma: Ayers et al., J Clin Invest
  2017;127(8):2930-2940. doi:10.1172/JCI91190.
    - 6-gene IFN-gamma set: VERIFIED — enumerated verbatim in the project
      FEATURE_CONTRACT_v2.md.
    - 18-gene list: canonical set as widely reproduced in the secondary
      literature; NOT verified against the primary Ayers supplement in this
      build (only the paper front-matter was retrievable — no body tables /
      gene enumeration). Use as the standard reconstruction, not a
      primary-verified original.
    - Exact GEP coefficients: UNVERIFIED (proprietary supplement not
      retrieved); equal-weight default used.
* T-effector & TGF-beta (pan-F-TBRS) gene sets: Mariathasan et al., Nature
  2018;554:544-548. doi:10.1038/nature25501 (IMvigor210). Gene sets as
  reproduced in that work's supplement / IMvigor210CoreBiologies package.

Real TPM source: this module runs on ANY gene x sample TPM matrix. In the
project pipeline the real matrix is `quant_gene_tpm.parquet` (Salmon TPM, per
HANDOFF_CONTRACT / FEATURE_CONTRACT_v2). This file ships a synthetic demo only.
"""
from __future__ import annotations

import logging
import warnings
from typing import Iterable, Literal, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gene signatures (HGNC symbols)
# ---------------------------------------------------------------------------

#: Ayers 2017 18-gene T-cell-inflamed GEP (the validated pembrolizumab
#: signature). NanoString-assay symbols. PROVENANCE: canonical gene list as
#: widely reproduced in the secondary literature; NOT verified against the
#: primary Ayers supplement in this session (only the paper front-matter was
#: retrievable, no body tables). Treat the list as the standard reconstruction,
#: not a primary-source-verified original.
GEP_TCELL_INFLAMED_GENES: tuple[str, ...] = (
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
)

#: Ayers 2017 6-gene IFN-gamma signature. VERIFIED against
#: FEATURE_CONTRACT_v2.md, which enumerates this 6-gene set verbatim.
IFNG_GENES: tuple[str, ...] = (
    "IDO1", "CXCL10", "CXCL9", "HLA-DRA", "IFNG", "STAT1",
)

#: Reconstructed T-cell-inflamed GEP weights (per-gene coefficients).
#: UNVERIFIED — reproduced from secondary sources, not the primary proprietary
#: supplement. Used only when weighting="published". Missing genes default to 0.
AYERS_GEP_WEIGHTS_PUBLISHED: dict[str, float] = {
    "CCL5": 0.008346, "CD27": 0.076189, "CD274": 0.042853, "CD276": 0.004313,
    "CD8A": 0.031021, "CMKLR1": 0.151253, "CXCL9": 0.074135, "CXCR6": 0.004313,
    "HLA-DQA1": 0.020091, "HLA-DRB1": 0.058806, "HLA-E": 0.058806,
    "IDO1": 0.060679, "LAG3": 0.123895, "NKG7": 0.075524, "PDCD1LG2": 0.003734,
    "PSMB10": 0.032999, "STAT1": 0.250229, "TIGIT": 0.084767,
}

#: Housekeeping genes for per-sample normalisation. The Ayers clinical assay
#: uses a proprietary 11-gene HK panel; the set below is the commonly-cited
#: reconstruction (stably-expressed reference genes). Any subset present in the
#: matrix is used to center each sample. Flagged UNVERIFIED against the assay.
AYERS_HOUSEKEEPING_GENES: tuple[str, ...] = (
    "STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "NRDE2",
    "UBB", "TBP", "SDHA", "ABCF1",
)

#: Mariathasan 2018 T-effector / IFN-gamma effector gene set (IMvigor210).
TEFF_GENES: tuple[str, ...] = (
    "CD8A", "EOMES", "PRF1", "IFNG", "CD274", "GZMA", "GZMB", "CXCL9",
    "CXCL10", "TBX21",
)

#: Mariathasan 2018 pan-F-TBRS (TGF-beta response) gene set (IMvigor210).
TGFB_GENES: tuple[str, ...] = (
    "ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCDC80",
    "COL5A1", "COL5A2", "COL6A3", "CTGF", "CTPS1", "RFLNB", "FBN1", "FN1",
    "FOXS1", "GREM1", "IGFBP3", "PMEPA1", "SGK1", "TGFBI", "TNS1", "TWIST1",
)

Harmonize = Literal["zscore", "rank", "none"]
Weighting = Literal["equal", "published"]


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def _to_log_tpm(tpm: pd.DataFrame) -> pd.DataFrame:
    """log2(TPM + 1). Rows = genes, cols = samples."""
    if (tpm.values < 0).any():
        raise ValueError("TPM matrix contains negative values; expected TPM >= 0.")
    return np.log2(tpm + 1.0)


def _resolve_genes(
    matrix_index: pd.Index, genes: Sequence[str], sig_name: str
) -> list[str]:
    """Intersect a signature gene list with the matrix, warn on missing."""
    present = [g for g in genes if g in matrix_index]
    missing = [g for g in genes if g not in matrix_index]
    if missing:
        logger.warning(
            "%s: %d/%d signature genes missing from matrix: %s",
            sig_name, len(missing), len(genes), ",".join(missing),
        )
    if not present:
        raise ValueError(
            f"{sig_name}: none of the {len(genes)} signature genes are present "
            "in the TPM matrix — cannot score."
        )
    return present


def _housekeeping_center(
    log_tpm: pd.DataFrame, housekeeping: Sequence[str] | None
) -> pd.DataFrame:
    """Subtract each SAMPLE's mean housekeeping log-expression (per-sample
    library-size / input-amount normalisation, as in the Ayers assay)."""
    if not housekeeping:
        return log_tpm
    hk_present = [g for g in housekeeping if g in log_tpm.index]
    if not hk_present:
        warnings.warn(
            "housekeeping normalisation requested but no housekeeping genes "
            "found in matrix; skipping HK-centering.",
            RuntimeWarning,
        )
        return log_tpm
    hk_mean = log_tpm.loc[hk_present].mean(axis=0)  # per-sample scalar
    return log_tpm.sub(hk_mean, axis=1)


def _sig_raw_value(
    log_tpm: pd.DataFrame,
    genes: Sequence[str],
    weights: dict[str, float] | None,
) -> pd.Series:
    """Weighted (or equal) mean of housekeeping-normalised log-TPM over the
    signature genes, per sample. Returns a per-sample Series."""
    block = log_tpm.loc[genes]  # genes x samples
    if weights is None:
        return block.mean(axis=0)
    w = pd.Series({g: weights.get(g, 0.0) for g in genes})
    if w.sum() == 0:
        raise ValueError("all published weights are zero for present genes.")
    return block.mul(w, axis=0).sum(axis=0) / w.sum()


def _harmonize(
    score: pd.Series, batches: pd.Series | None, method: Harmonize
) -> pd.Series:
    """Standardise a score WITHIN each batch (platform robustness)."""
    if method == "none":
        return score
    if batches is None:
        batches = pd.Series("_all_", index=score.index)

    def _one(s: pd.Series) -> pd.Series:
        if method == "zscore":
            sd = s.std(ddof=0)
            if sd == 0 or np.isnan(sd):
                return pd.Series(0.0, index=s.index)
            return (s - s.mean()) / sd
        elif method == "rank":
            # rank -> (0,1) percentile, then to N(0,1)-ish centered [-0.5,0.5]
            r = s.rank(method="average")
            return (r - 0.5) / len(s) - 0.5
        raise ValueError(f"unknown harmonize method {method!r}")

    return score.groupby(batches, group_keys=False).apply(_one).reindex(score.index)


# ---------------------------------------------------------------------------
# Public scorers
# ---------------------------------------------------------------------------

def _prep(
    tpm: pd.DataFrame, housekeeping: Sequence[str] | None
) -> pd.DataFrame:
    log_tpm = _to_log_tpm(tpm)
    return _housekeeping_center(log_tpm, housekeeping)


def score_gep_tcell_inflamed(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    weighting: Weighting = "equal",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """Ayers 2017 18-gene T-cell-inflamed GEP score. -> gep_tcell_inflamed.

    Parameters
    ----------
    tpm : genes (index) x samples (columns) TPM matrix.
    batches : per-sample batch/cohort labels (index == tpm.columns). Score is
        standardised within each batch when harmonize != "none".
    harmonize : "zscore" (default), "rank", or "none".
    weighting : "equal" (default, reproducible) or "published" (reconstructed
        Ayers coefficients — UNVERIFIED).
    housekeeping : housekeeping genes for per-sample centering; None to skip.
    """
    log_tpm = _prep(tpm, housekeeping)
    genes = _resolve_genes(log_tpm.index, GEP_TCELL_INFLAMED_GENES, "gep_tcell_inflamed")
    weights = None
    if weighting == "published":
        weights = {g: AYERS_GEP_WEIGHTS_PUBLISHED[g] for g in genes}
    raw = _sig_raw_value(log_tpm, genes, weights)
    return _harmonize(raw, batches, harmonize).rename("gep_tcell_inflamed")


def score_ifng(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """Ayers 2017 6-gene IFN-gamma signature score (equal-weight mean, as
    published). -> ifng_score."""
    log_tpm = _prep(tpm, housekeeping)
    genes = _resolve_genes(log_tpm.index, IFNG_GENES, "ifng_score")
    raw = _sig_raw_value(log_tpm, genes, weights=None)
    return _harmonize(raw, batches, harmonize).rename("ifng_score")


def score_teff_tgfb_balance(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """Mariathasan 2018 T-effector minus TGF-beta (pan-F-TBRS) balance.
    -> teff_tgfb_balance.

    High = T-effector-dominant (immune-active); low = TGF-beta/stromal-dominant
    (associated with ICB resistance). Each arm is an equal-weight mean of
    housekeeping-normalised log-TPM; the two arms are each standardised within
    batch BEFORE subtracting, so the balance is on a common per-batch scale.
    """
    log_tpm = _prep(tpm, housekeeping)
    teff_genes = _resolve_genes(log_tpm.index, TEFF_GENES, "teff (Mariathasan)")
    tgfb_genes = _resolve_genes(log_tpm.index, TGFB_GENES, "tgfb pan-F-TBRS (Mariathasan)")
    teff_raw = _sig_raw_value(log_tpm, teff_genes, weights=None)
    tgfb_raw = _sig_raw_value(log_tpm, tgfb_genes, weights=None)
    # standardise each arm within batch first, then take the difference so the
    # balance is not dominated by whichever arm has the larger dynamic range.
    teff_z = _harmonize(teff_raw, batches, harmonize if harmonize != "none" else "zscore")
    tgfb_z = _harmonize(tgfb_raw, batches, harmonize if harmonize != "none" else "zscore")
    return (teff_z - tgfb_z).rename("teff_tgfb_balance")


# ---------------------------------------------------------------------------
# Batch entrypoint -> feature table
# ---------------------------------------------------------------------------

def score_all(
    tpm: pd.DataFrame,
    sample_meta: pd.DataFrame,
    *,
    harmonize: Harmonize = "zscore",
    weighting: Weighting = "equal",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
    batch_col: str = "cohort",
) -> pd.DataFrame:
    """Compute all three baseline GEP features for every sample.

    Parameters
    ----------
    tpm : genes x samples TPM matrix. Columns must match sample_meta index or
        its run_accession.
    sample_meta : per-sample metadata; MUST contain 'run_accession' and the
        batch column (default 'cohort'). Its row order defines output order.
        Indexed by the same sample ids as tpm.columns.
    batch_col : column in sample_meta giving the within-batch standardisation
        unit (default 'cohort').

    Returns
    -------
    DataFrame keyed on run_accession + cohort with columns
    [gep_tcell_inflamed, ifng_score, teff_tgfb_balance].
    """
    if "run_accession" not in sample_meta.columns:
        raise ValueError("sample_meta must contain a 'run_accession' column.")
    if batch_col not in sample_meta.columns:
        raise ValueError(f"sample_meta must contain the batch column {batch_col!r}.")

    # align columns of tpm to sample_meta
    samples = list(sample_meta.index)
    missing_cols = [s for s in samples if s not in tpm.columns]
    if missing_cols:
        raise ValueError(
            f"{len(missing_cols)} samples in sample_meta absent from TPM matrix: "
            f"{missing_cols[:5]}..."
        )
    tpm = tpm[samples]
    batches = sample_meta[batch_col]

    gep = score_gep_tcell_inflamed(
        tpm, batches, harmonize=harmonize, weighting=weighting, housekeeping=housekeeping
    )
    ifng = score_ifng(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)
    teff = score_teff_tgfb_balance(
        tpm, batches, harmonize=harmonize, housekeeping=housekeeping
    )

    out = pd.DataFrame(
        {
            "run_accession": sample_meta["run_accession"].values,
            "cohort": sample_meta[batch_col].values,
            "gep_tcell_inflamed": gep.values,
            "ifng_score": ifng.values,
            "teff_tgfb_balance": teff.values,
        },
        index=sample_meta.index,
    )
    return out.reset_index(drop=True)


FEATURE_COLUMNS = ("gep_tcell_inflamed", "ifng_score", "teff_tgfb_balance")
