"""
src/features.py — consume the de-novo phenotype matrices produced by the
pipeline session, per docs/HANDOFF_CONTRACT.md.

The modeling layer never runs STAR/nf-core. It reads per-sample feature
matrices from ``results/features/`` (rows = ENA run_accession, cols = numeric
features, first two cols run_accession + cohort) and joins them to the frozen
clinical/response frame from ``src/data.py``.

Degrades gracefully: ``load_available()`` returns only the matrices that exist,
so the modeling code can run the moment the first phenotype lands and expand as
more arrive.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FEAT_DIR = REPO / "results" / "features"
RESULTS = REPO / "results"

# phenotype id -> expected filename (see HANDOFF_CONTRACT.md)
EXPECTED = {
    "quant_gene_tpm": "quant_gene_tpm.parquet",
    "quant_tx_tpm": "quant_tx_tpm.parquet",
    "splicing_psi": "splicing_psi.parquet",
    "intron_retention": "intron_retention.parquet",
    "rna_editing_aei": "rna_editing_aei.parquet",
    "te_locus": "te_locus.parquet",
    "te_family": "te_family.parquet",
    "fusion_burden": "fusion_burden.parquet",
}
KEY = "run_accession"
META_COLS = [KEY, "cohort"]


def _read(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def load_available(feat_dir: Path = FEAT_DIR) -> dict[str, pd.DataFrame]:
    """Return {phenotype: DataFrame} for every matrix currently present.

    Accepts .parquet or .csv. Validates the sample-key column exists; warns
    (does not raise) on a matrix that lacks it so a malformed drop is visible.
    """
    out: dict[str, pd.DataFrame] = {}
    if not feat_dir.exists():
        return out
    for pheno, fname in EXPECTED.items():
        for cand in (feat_dir / fname,
                     feat_dir / fname.replace(".parquet", ".csv")):
            if cand.exists():
                df = _read(cand)
                if KEY not in df.columns:
                    print(f"[features] WARNING {cand.name}: no '{KEY}' column; skipped")
                else:
                    out[pheno] = df
                break
    return out


def status(feat_dir: Path = FEAT_DIR) -> pd.DataFrame:
    """One row per expected phenotype: present?, n_samples, n_features."""
    avail = load_available(feat_dir)
    rows = []
    for pheno, fname in EXPECTED.items():
        if pheno in avail:
            df = avail[pheno]
            nfeat = df.shape[1] - sum(c in df.columns for c in META_COLS)
            rows.append({"phenotype": pheno, "present": True,
                         "n_samples": len(df), "n_features": nfeat})
        else:
            rows.append({"phenotype": pheno, "present": False,
                         "n_samples": 0, "n_features": 0})
    return pd.DataFrame(rows)


def join_clinical(avail: dict[str, pd.DataFrame],
                  clinical_key: str = "run_accession") -> pd.DataFrame:
    """Join available feature matrices to clinical/response labels.

    Uses data/manifests/selection_manifest.csv (has run_accession + response)
    as the clinical bridge, since the frozen analysis_frame is keyed by
    cBioPortal sampleId, not ENA accession.
    """
    manifest = pd.read_csv(REPO / "data" / "manifests" / "selection_manifest.csv")
    keep = [c for c in ["cohort", "run_accession", "patient_id", "timepoint",
                         "recist", "resp_NR", "therapy_clin", "os_days", "vital"]
            if c in manifest.columns]
    base = manifest[keep].copy()
    for pheno, df in avail.items():
        feat_cols = [c for c in df.columns if c not in META_COLS]
        renamed = {c: f"{pheno}::{c}" for c in feat_cols}
        base = base.merge(df.rename(columns=renamed)[[KEY] + list(renamed.values())],
                          on=KEY, how="left")
    return base


if __name__ == "__main__":
    print(status().to_string(index=False))
    avail = load_available()
    print("\navailable phenotypes:", sorted(avail))
