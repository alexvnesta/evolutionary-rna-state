"""
src/data.py — build the harmonized per-sample analysis frame.

Joins the five cBioPortal iAtlas-harmonized clinical pulls
(``data/cbioportal/*.clinical.tsv``) into one tidy table with a stable schema:

  * phenotype block  — mechanism-resolved neoantigen loads that proxy the
                       thesis RNA-phenotypes (SPLICE / ERV / FUSION), plus the
                       mutational-neoantigen categories (INDEL / SNV) used as
                       burden controls.
  * confounders      — TMB / mutation count (burden) and, for DFCI only,
                       purity / ploidy / clonal-mutation counts.
  * response         — RESPONDER (bool) and RECIST RESPONSE.
  * survival         — OS / PFS months + status.
  * benchmark        — TIDE_RESPONDER (published response-signature call).

Authored in-session (hackathon clean-slate). Metadata-only inputs; no raw reads.

Run as a script to write ``results/analysis_frame.parquet`` and a coverage/
missingness figure ``results/fig_coverage.png``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[1]
CBIO_DIR = REPO / "data" / "cbioportal"
RESULTS = REPO / "results"

# cohort id -> (filename glob stem, cBioPortal study id)
COHORTS = {
    "gide2019": "mel_iatlas_gide_2019",
    "riaz2017": "mel_iatlas_riaz_nivolumab_2017",
    "hugo2016": "mel_iatlas_hugo_ucla_2016",
    "liu2019": "mel_iatlas_liu_2019",
    "dfci2019": "mel_dfci_2019",
}

# ---------------------------------------------------------------------------
# Column schema (source column names as pulled from cBioPortal)
# ---------------------------------------------------------------------------
# Mechanism-resolved neoantigen loads. SPLICE / ERV / FUSION are the three
# thesis RNA-phenotypes with a per-sample proxy; INDEL / SNV are mutational
# categories retained as burden controls (not RNA-phenotypes themselves).
PHENO_RNA = ["SPLICE_NEOANTIGEN", "ERV_NEOANTIGEN", "FUSION_NEOANTIGEN"]
PHENO_MUT = ["INDEL_NEOANTIGEN", "SNV_NEOANTIGEN"]
PHENO_ALL = PHENO_RNA + PHENO_MUT

# Burden confounders ("the field standard to beat").
BURDEN = ["TMB_NONSYNONYMOUS", "MUTATION_COUNT"]

# Composition / genome-instability covariates (DFCI carries these).
INSTABILITY = ["PURITY", "PLOIDY", "MUTS_CLONAL", "MUTS_SUBCLONAL",
               "HETEROGENEITY", "CNA_PROP"]

# Response, survival, benchmark.
RESPONSE = ["RESPONDER", "RESPONSE", "CLINICAL_BENEFIT"]
SURVIVAL = ["OS_MONTHS", "OS_STATUS", "PFS_MONTHS", "PFS_STATUS"]
BENCH = ["TIDE_RESPONDER"]

META = ["patientId", "SAMPLE_TREATMENT", "CANCER_TYPE_DETAILED",
        "ICI_TARGET", "ICI_RX", "SEX", "AGE_AT_DIAGNOSIS"]

KEEP = (["sampleId", "cohort"] + META + PHENO_ALL + BURDEN + INSTABILITY
        + RESPONSE + SURVIVAL + BENCH)


# ---------------------------------------------------------------------------
# Loading / harmonization
# ---------------------------------------------------------------------------
def _cohort_path(cohort: str) -> Path:
    hits = sorted(CBIO_DIR.glob(f"{cohort}_*.clinical.tsv"))
    if not hits:
        raise FileNotFoundError(f"no clinical tsv for {cohort} in {CBIO_DIR}")
    return hits[0]


def load_cohort(cohort: str) -> pd.DataFrame:
    """Load one cohort, add any missing schema columns as NA, tag cohort."""
    df = pd.read_csv(_cohort_path(cohort), sep="\t")
    df["cohort"] = cohort
    # SAMPLE_TREATMENT absent (DFCI) => treat all rows as pretreatment.
    if "SAMPLE_TREATMENT" not in df.columns:
        df["SAMPLE_TREATMENT"] = "Pre"
    for col in KEEP:
        if col not in df.columns:
            df[col] = np.nan
    return df[KEEP].copy()


def load_all() -> pd.DataFrame:
    """Concatenate all cohorts into one long analysis frame."""
    frames = [load_cohort(c) for c in COHORTS]
    df = pd.concat(frames, ignore_index=True)
    # Normalize RESPONDER to nullable boolean.
    df["RESPONDER"] = df["RESPONDER"].map(
        {True: True, False: False, "True": True, "False": False,
         "Yes": True, "No": False}
    ).astype("boolean")
    return df


def freeze_pretreatment(df: pd.DataFrame) -> pd.DataFrame:
    """Freeze the pretreatment analysis set.

    Keep pretreatment ('Pre') samples only. When a patient has multiple
    pretreatment samples, keep the first (stable order) so the unit of
    analysis is one row per patient-baseline.
    """
    pre = df[df["SAMPLE_TREATMENT"].astype(str).str.lower().eq("pre")].copy()
    pre = pre.sort_values(["cohort", "patientId", "sampleId"])
    pre = pre.drop_duplicates(subset=["cohort", "patientId"], keep="first")
    return pre.reset_index(drop=True)


def phenotype_cohorts(pre: pd.DataFrame) -> pd.DataFrame:
    """Subset to cohorts that actually carry the RNA-phenotype proxy loads.

    DFCI has no neoantigen-mechanism columns (covariate-only cohort), so it is
    excluded from the phenotype analysis set and used separately for the
    purity / clonal-instability covariate context.
    """
    has_pheno = pre[PHENO_RNA].notna().any(axis=1)
    return pre[has_pheno].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Coverage / missingness report
# ---------------------------------------------------------------------------
def coverage_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Per-cohort non-null fraction for every analytic field."""
    fields = PHENO_ALL + BURDEN + INSTABILITY + ["RESPONDER", "RESPONSE",
             "OS_MONTHS", "PFS_MONTHS", "TIDE_RESPONDER"]
    rows = {}
    for cohort, sub in df.groupby("cohort"):
        rows[cohort] = {f: sub[f].notna().mean() for f in fields}
    cov = pd.DataFrame(rows).T
    cov["n"] = df.groupby("cohort").size()
    return cov[["n"] + fields]


def build(write: bool = True) -> pd.DataFrame:
    """End-to-end: load, freeze pretreatment, tag phenotype set, persist."""
    RESULTS.mkdir(exist_ok=True)
    full = load_all()
    pre = freeze_pretreatment(full)
    pheno = phenotype_cohorts(pre)
    pre["in_phenotype_set"] = pre["sampleId"].isin(pheno["sampleId"])
    if write:
        pre.to_parquet(RESULTS / "analysis_frame.parquet", index=False)
        cov = coverage_matrix(pre)
        cov.to_csv(RESULTS / "coverage_matrix.csv")
        meta = {
            "n_total_samples": int(len(full)),
            "n_pretreatment": int(len(pre)),
            "n_phenotype_set": int(len(pheno)),
            "phenotype_cohorts": sorted(pheno["cohort"].unique().tolist()),
            "phenotype_rna_cols": PHENO_RNA,
            "phenotype_mut_cols": PHENO_MUT,
            "burden_cols": BURDEN,
            "per_cohort_pretreatment_n":
                pre.groupby("cohort").size().to_dict(),
            "per_cohort_phenotype_n":
                pheno.groupby("cohort").size().to_dict(),
        }
        (RESULTS / "analysis_set_meta.json").write_text(json.dumps(meta, indent=2))
    return pre


if __name__ == "__main__":
    pre = build(write=True)
    print(pre.groupby("cohort").size())
    print("phenotype set n =", int(pre["in_phenotype_set"].sum()))
