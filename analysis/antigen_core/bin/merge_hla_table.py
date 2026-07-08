#!/usr/bin/env python
"""
bin/merge_hla_table.py — merge per-sample arcasHLA genotype JSONs into the
cohort tidy HLA table (contract format). Called by the MERGE_HLA_TABLE
Nextflow process.

Keys each genotype JSON (basename = <run_accession>.genotype.json) to its
run_accession + cohort via --sample-map (a CSV with run_accession,cohort),
then applies the shared heterozygosity logic from
analysis/antigen_core/hla_typing.py.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# make the antigen_core package importable whether run from repo or Nextflow work dir
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

from hla_typing import parse_arcashla_genotype, summarize_genotype, build_hla_table  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotypes", nargs="+", required=True,
                    help="per-sample *.genotype.json files")
    ap.add_argument("--sample-map", required=True,
                    help="CSV with columns run_accession,cohort")
    ap.add_argument("--tool", default="arcasHLA")
    ap.add_argument("--tool-version", default="")
    ap.add_argument("--out", required=True, help="output parquet path")
    args = ap.parse_args()

    smap = pd.read_csv(args.sample_map)
    cohort_by_run = dict(zip(smap["run_accession"], smap["cohort"]))

    rows = []
    for gj in args.genotypes:
        run = Path(gj).name.split(".genotype.json")[0]
        cohort = cohort_by_run.get(run, "")
        genotype = parse_arcashla_genotype(gj)
        rows.append(summarize_genotype(
            genotype, run, cohort,
            tool=args.tool, tool_version=args.tool_version,
        ))

    df = build_hla_table(rows)
    if args.out.endswith(".parquet"):
        df.to_parquet(args.out, index=False)
    else:
        df.to_csv(args.out, index=False)
    print(f"[merge_hla_table] wrote {len(df)} samples -> {args.out}")


if __name__ == "__main__":
    main()
