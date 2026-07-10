"""
analysis/baseline/apm_scores.py
===============================

BASELINE BUCKET — antigen-presentation-machinery (APM) / HLA expression score.

A named, batch-robust gene-expression score of the MHC class-I antigen
processing-and-presentation machinery, computed from a gene x sample TPM
matrix using EXACTLY the same conventions as gep_scores.py (log2(TPM+1),
per-sample housekeeping centering, within-cohort z-score harmonisation). This
is a *proven baseline* axis: loss / low expression of the class-I APM is one of
the most reproducible pre-treatment resistance markers to checkpoint blockade,
and it is mechanistically the gate on whether any neoantigen — canonical or
RNA-derived — can be presented at all.

Three named scores are exposed:

    apm_class1        -- MHC class-I processing+presentation machinery
    apm_class2        -- MHC class-II machinery (companion; overlaps immune
                         infiltrate more than class-I — see CAVEAT)
    b2m_hla_abc       -- the minimal HLA-A/B/C + B2M "presentation floor",
                         the sub-score most specific to tumour-cell-intrinsic
                         presentation competence

All are floats, one value per sample, keyed on (run_accession, cohort).

Why a SEPARATE class-I score rather than folding into GEP
---------------------------------------------------------
The Ayers T-cell-inflamed GEP already contains HLA-E and PSMB10 (and class-II
HLA-DQA1/HLA-DRB1). The class-I APM score here is deliberately built from the
tumour-cell-intrinsic processing/presentation components with MINIMAL GEP
overlap (only HLA-E and PSMB10 are shared, and both can be dropped via
`exclude_gep_overlap=True`) so that, in the downstream evaluation, "APM adds
over GEP" is a question about presentation competence and not a re-statement of
the same 18 genes. This is the proxy-circularity guard applied at design time.

CAVEAT (composition confounding)
--------------------------------
HLA class-I and especially class-II expression in a *bulk* tumour sample is a
mixture of tumour-cell-intrinsic expression and the APM carried by infiltrating
immune cells (which are HLA-high). A high bulk APM score can therefore reflect
"more immune cells" rather than "tumour cells competent to present." The
class-II arm is the more infiltrate-driven of the two. Any association with
response MUST be checked for composition confounding (signature-rigour-harness)
before interpretation — the same discipline the project applies to every
signature.

Gene-set provenance
-------------------
There is no single canonical primary-source "APM score gene list"; multiple ICB
papers use overlapping component sets assembled from the antigen-presentation
pathway. The lists below are CURATED from the well-established components of the
MHC-I / MHC-II antigen processing-and-presentation pathway (HLA loci, peptide
transporters, the immunoproteasome, ER peptide-loading complex, ER
aminopeptidases, and the CIITA/NLRC5 master transactivators). Individual gene
membership in the pathway is textbook / HGNC-verifiable, but this exact
composite is a RECONSTRUCTION, not a verbatim copy of one published signature.
Treat it as the standard component set, flagged accordingly — the same posture
gep_scores.py takes for the 18-gene GEP list.

This module intentionally reuses gep_scores' internal helpers so there is ONE
implementation of log-transform / housekeeping-centering / harmonisation.
"""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from analysis.baseline import gep_scores as _gs
from analysis.baseline.gep_scores import Harmonize, AYERS_HOUSEKEEPING_GENES

# ---------------------------------------------------------------------------
# Gene signatures (HGNC symbols)
# ---------------------------------------------------------------------------

#: MHC class-I antigen processing & presentation machinery. CURATED pathway
#: component set (see module provenance). Components:
#:   - HLA class-I heavy chains + invariant light chain: HLA-A/B/C, B2M
#:   - peptide transporter into the ER: TAP1, TAP2, TAPBP (tapasin)
#:   - ER peptide-loading complex chaperones: CALR, CANX, PDIA3
#:   - immunoproteasome catalytic subunits: PSMB8, PSMB9, PSMB10
#:   - ER aminopeptidases (peptide trimming): ERAP1, ERAP2
#:   - master transactivator of class-I: NLRC5
APM_CLASS1_GENES: tuple[str, ...] = (
    "HLA-A", "HLA-B", "HLA-C", "B2M",
    "TAP1", "TAP2", "TAPBP",
    "CALR", "CANX", "PDIA3",
    "PSMB8", "PSMB9", "PSMB10",
    "ERAP1", "ERAP2",
    "NLRC5",
)

#: Genes in APM_CLASS1 that also appear in the Ayers T-cell-inflamed GEP.
#: Dropped when exclude_gep_overlap=True so the class-I score is orthogonal to
#: the GEP floor it will be compared against.
APM_CLASS1_GEP_OVERLAP: tuple[str, ...] = ("PSMB10",)  # HLA-E is in GEP but not in this class-I set

#: MHC class-II machinery (companion score). CURATED pathway component set:
#:   - classical class-II heterodimers: HLA-DRA/DRB1, HLA-DPA1/DPB1, HLA-DQA1/DQB1
#:   - class-II peptide-editing molecules: HLA-DMA, HLA-DMB
#:   - invariant chain: CD74
#:   - master transactivator of class-II: CIITA
APM_CLASS2_GENES: tuple[str, ...] = (
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "HLA-DMA", "HLA-DMB", "CD74", "CIITA",
)

#: Minimal "presentation floor": the HLA-A/B/C heavy chains + B2M. This is the
#: sub-score most specific to tumour-cell-intrinsic class-I presentation
#: competence (B2M loss / HLA-ABC loss is the textbook hard ICB-resistance
#: lesion). No GEP overlap.
B2M_HLA_ABC_GENES: tuple[str, ...] = ("HLA-A", "HLA-B", "HLA-C", "B2M")


# ---------------------------------------------------------------------------
# Public scorers (thin wrappers over the shared gep_scores machinery)
# ---------------------------------------------------------------------------

def _score_set(
    tpm: pd.DataFrame,
    genes: Sequence[str],
    name: str,
    batches: pd.Series | None,
    harmonize: Harmonize,
    housekeeping: Sequence[str] | None,
) -> pd.Series:
    """Equal-weight mean of housekeeping-normalised log-TPM over `genes`,
    harmonised within batch — identical pipeline to gep_scores' scorers."""
    log_tpm = _gs._prep(tpm, housekeeping)
    present = _gs._resolve_genes(log_tpm.index, genes, name)
    raw = _gs._sig_raw_value(log_tpm, present, weights=None)
    return _gs._harmonize(raw, batches, harmonize).rename(name)


def score_apm_class1(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
    exclude_gep_overlap: bool = False,
) -> pd.Series:
    """MHC class-I APM expression score. -> apm_class1.

    exclude_gep_overlap: drop the genes shared with the Ayers GEP so the score
    is orthogonal to the GEP floor (default False; set True for the
    proxy-circularity-guarded comparison).
    """
    genes = APM_CLASS1_GENES
    if exclude_gep_overlap:
        genes = tuple(g for g in genes if g not in set(APM_CLASS1_GEP_OVERLAP))
    return _score_set(tpm, genes, "apm_class1", batches, harmonize, housekeeping)


def score_apm_class2(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """MHC class-II APM expression score (companion). -> apm_class2."""
    return _score_set(tpm, APM_CLASS2_GENES, "apm_class2", batches, harmonize, housekeeping)


def score_b2m_hla_abc(
    tpm: pd.DataFrame,
    batches: pd.Series | None = None,
    *,
    harmonize: Harmonize = "zscore",
    housekeeping: Sequence[str] | None = AYERS_HOUSEKEEPING_GENES,
) -> pd.Series:
    """Minimal HLA-A/B/C + B2M presentation-floor score. -> b2m_hla_abc."""
    return _score_set(tpm, B2M_HLA_ABC_GENES, "b2m_hla_abc", batches, harmonize, housekeeping)


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
    exclude_gep_overlap: bool = False,
) -> pd.DataFrame:
    """Compute all APM features for every sample.

    tpm : genes(symbol) x samples matrix. sample_meta : per-sample metadata
    with 'run_accession' and the batch column; its row order defines output
    order. Mirrors gep_scores.score_all's contract exactly.
    """
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

    c1 = score_apm_class1(tpm, batches, harmonize=harmonize, housekeeping=housekeeping,
                          exclude_gep_overlap=exclude_gep_overlap)
    c2 = score_apm_class2(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)
    floor = score_b2m_hla_abc(tpm, batches, harmonize=harmonize, housekeeping=housekeeping)

    out = pd.DataFrame(
        {
            "run_accession": sample_meta["run_accession"].values,
            "cohort": sample_meta[batch_col].values,
            "apm_class1": c1.values,
            "apm_class2": c2.values,
            "b2m_hla_abc": floor.values,
        },
        index=sample_meta.index,
    )
    return out.reset_index(drop=True)


FEATURE_COLUMNS = ("apm_class1", "apm_class2", "b2m_hla_abc")
