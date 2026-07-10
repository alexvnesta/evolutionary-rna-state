"""
analysis/baseline/tls_bcell_scores.py
======================================

BASELINE BUCKET — tertiary-lymphoid-structure (TLS) / B-cell signature score.

Named, batch-robust gene-expression scores of the intratumoural B-cell / TLS
compartment, computed from a gene x sample TPM matrix with EXACTLY the
gep_scores.py conventions (log2(TPM+1), per-sample housekeeping centering,
within-cohort z-score harmonisation).

Motivation: tumour-infiltrating B cells and tertiary lymphoid structures are
reproducible PRE-treatment predictors of checkpoint-blockade response that are
largely ORTHOGONAL to the T-cell-inflamed axis (melanoma, sarcoma, RCC). This
is also the population-level, RNA-observable reading of B-cell somatic
hypermutation / affinity maturation — the same biology the TRUST4 BCR module
measures at the repertoire level.

Three named scores:

    bcell_lineage     -- core B-cell / plasma-cell lineage abundance
    tls_chemokine     -- the TLS-associated chemokine programme
    tls_imprint       -- combined B-cell + TLS-chemokine imprint (mean of the
                         two standardised arms)

All are floats, one value per sample, keyed on (run_accession, cohort).

CAVEAT (composition confounding)
--------------------------------
These are, by construction, ABUNDANCE readouts of a cell compartment: a high
score largely means "more B cells / more TLS-like organisation in the bulk
sample." That is exactly what makes them predictive, but it also means any
association with response is a composition signal and must be evaluated as such
(signature-rigour-harness) — the orthogonality to interrogate is B-cell/TLS vs
the T-cell GEP, not vs immune infiltrate in general.

Gene-set provenance
-------------------
* bcell_lineage: canonical B-cell and plasma-cell LINEAGE markers (CD20/MS4A1,
  CD19, CD79A/B, the TACI receptor TNFRSF13B, and the plasma-cell markers MZB1,
  DERL3, TNFRSF17/BCMA, SDC1/CD138). Individual gene identity is textbook /
  HGNC-verifiable. Assembled here as a lineage panel; NOT a verbatim copy of a
  single published "B-cell signature."
* tls_chemokine: the widely-used TLS-associated 12-chemokine programme
  (Coppola/Messina-type). Membership below is the commonly-reproduced set and
  is flagged a RECONSTRUCTION — it was NOT verified against the primary
  supplement in this build. Use as the standard reconstruction, and verify
  against the primary source before any membership claim enters a writeup
  (same posture gep_scores.py takes for the 18-gene GEP list).

This module reuses gep_scores' internal helpers so there is ONE implementation
of the transform / centering / harmonisation pipeline.
"""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from analysis.baseline import gep_scores as _gs
from analysis.baseline.gep_scores import Harmonize, AYERS_HOUSEKEEPING_GENES

# ---------------------------------------------------------------------------
# Gene signatures (HGNC symbols)
# ---------------------------------------------------------------------------

#: Canonical B-cell + plasma-cell lineage markers (see provenance).
BCELL_LINEAGE_GENES: tuple[str, ...] = (
    "MS4A1", "CD19", "CD79A", "CD79B", "TNFRSF13B",   # B cell
    "MZB1", "DERL3", "TNFRSF17", "SDC1",              # plasma cell
)

#: TLS-associated chemokine programme (12-chemokine, Coppola/Messina-type).
#: RECONSTRUCTION — commonly-reproduced membership, not primary-verified here.
TLS_CHEMOKINE_GENES: tuple[str, ...] = (
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
    "CXCL9", "CXCL10", "CXCL11", "CXCL13",
)

#: Genes shared with the Ayers T-cell-inflamed GEP (CCL5, CXCL9). Reported so
#: the B-cell/TLS-vs-GEP orthogonality can be interpreted honestly; the TLS
#: chemokine programme legitimately includes T-cell-recruiting chemokines, so
#: these are NOT dropped by default.
TLS_CHEMOKINE_GEP_OVERLAP: tuple[str, ...] = ("CCL5", "CXCL9")


def _score_set(
    tpm: pd.DataFrame,
    genes: Sequence[str],
    name: str,
    batches: pd.Series | None,
    harmonize: Harmonize,
    housekeeping: Sequence[str] | None,
) -> pd.Series:
    log_tpm = _gs._prep(tpm, housekeeping)
    present = _gs._resolve_genes(log_tpm.index, genes, name)
    raw = _gs._sig_raw_value(log_tpm, present, weights=None)
    return _gs._harmonize(raw, batches, harmonize).rename(name)


def score_bcell_lineage(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """B-cell / plasma-cell lineage abundance score. -> bcell_lineage."""
    return _score_set(tpm, BCELL_LINEAGE_GENES, "bcell_lineage", batches, harmonize, housekeeping)


def score_tls_chemokine(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """TLS-associated 12-chemokine programme score. -> tls_chemokine."""
    return _score_set(tpm, TLS_CHEMOKINE_GENES, "tls_chemokine", batches, harmonize, housekeeping)


def score_tls_imprint(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """Combined B-cell + TLS-chemokine imprint: mean of the two within-batch
    standardised arms (so neither arm's dynamic range dominates). -> tls_imprint.
    """
    b = score_bcell_lineage(tpm, batches, harmonize=harmonize if harmonize != "none" else "zscore",
                            housekeeping=housekeeping)
    c = score_tls_chemokine(tpm, batches, harmonize=harmonize if harmonize != "none" else "zscore",
                            housekeeping=housekeeping)
    return ((b + c) / 2.0).rename("tls_imprint")


# ---------------------------------------------------------------------------
# Batch entrypoint -> feature table
# ---------------------------------------------------------------------------

def score_all(
    tpm: pd.DataFrame,
    sample_meta: pd.DataFrame,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
    batch_col: str = "cohort",
) -> pd.DataFrame:
    """Compute all TLS/B-cell features for every sample. Contract mirrors
    gep_scores.score_all / apm_scores.score_all exactly."""
    if "run_accession" not in sample_meta.columns:
        raise ValueError("sample_meta must contain a 'run_accession' column.")
    if batch_col not in sample_meta.columns:
        raise ValueError(f"sample_meta must contain the batch column {batch_col!r}.")

    samples = list(sample_meta.index)
    missing_cols = [s for s in samples if s not in tpm.columns]
    if missing_cols:
        raise ValueError(
            f"{len(missing_cols)} samples in sample_meta absent from TPM matrix: "
            f"{missing_cols[:5]}..."
        )
    tpm = tpm[samples]
    batches = sample_meta[batch_col]

    bl = score_bcell_lineage(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)
    tc = score_tls_chemokine(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)
    ti = score_tls_imprint(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)

    out = pd.DataFrame(
        {
            "run_accession": sample_meta["run_accession"].values,
            "cohort": sample_meta[batch_col].values,
            "bcell_lineage": bl.values,
            "tls_chemokine": tc.values,
            "tls_imprint": ti.values,
        },
        index=sample_meta.index,
    )
    return out.reset_index(drop=True)


FEATURE_COLUMNS = ("bcell_lineage", "tls_chemokine", "tls_imprint")
