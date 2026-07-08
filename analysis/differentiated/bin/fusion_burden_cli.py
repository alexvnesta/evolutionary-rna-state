#!/usr/bin/env python
"""
bin/fusion_burden_cli.py — per-sample fusion-neoantigen-burden driver, called by
the FUSION_NEOANTIGEN_BURDEN Nextflow process.

Reads ONE caller's fusion TSV (Arriba XOR STAR-Fusion — the fixed-caller
requirement), pulls the sample's 6 HLA-I alleles from the cohort HLA table by
run_accession, and writes a single contract-shaped feature-row CSV via
analysis/differentiated/fusion_antigen.py.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# make the differentiated + antigen_core packages importable whether run from
# the repo or a Nextflow work dir (mirrors bin/merge_hla_table.py)
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))   # analysis/differentiated
sys.path.insert(0, str(_HERE.parents[2]))   # analysis  (-> analysis.antigen_core)

import pandas as pd  # noqa: E402

import fusion_antigen as fa  # noqa: E402

ALLELE_COLS = ["HLA_A_1", "HLA_A_2", "HLA_B_1", "HLA_B_2", "HLA_C_1", "HLA_C_2"]


def hla_for_run(hla_table: str, run_accession: str) -> list[str]:
    df = pd.read_parquet(hla_table) if str(hla_table).endswith(".parquet") \
        else pd.read_csv(hla_table)
    hit = df[df["run_accession"] == run_accession]
    if hit.empty:
        raise SystemExit(f"run_accession {run_accession} not in HLA table {hla_table}")
    row = hit.iloc[0]
    return [str(row[c]) for c in ALLELE_COLS if c in row and pd.notna(row[c])]


def main() -> None:
    ap = argparse.ArgumentParser()
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--arriba", help="Arriba fusions.tsv")
    grp.add_argument("--starfusion", help="STAR-Fusion coding-effect TSV")
    ap.add_argument("--hla-table", required=True, help="hla_typing.parquet/csv")
    ap.add_argument("--run-accession", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--caller-version", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    hla = hla_for_run(args.hla_table, args.run_accession)
    row = fa.fusion_features_for_sample(
        run_accession=args.run_accession,
        cohort=args.cohort,
        hla_alleles=hla,
        arriba_tsv=args.arriba,
        starfusion_tsv=args.starfusion,
        caller_version=args.caller_version,
    )
    fa.build_fusion_feature_table([row]).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
