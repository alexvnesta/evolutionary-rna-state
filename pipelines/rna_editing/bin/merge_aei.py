#!/usr/bin/env python3
"""
merge_aei.py -- concatenate per-sample AEI TSVs into one cohort table.
Usage: merge_aei.py --out cohort_aei.tsv S1.aei.tsv S2.aei.tsv ...
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("tsvs", nargs="+")
    args = ap.parse_args()
    df = pd.concat([pd.read_csv(t, sep="\t") for t in args.tsvs], ignore_index=True)
    df = df.sort_values("sample").reset_index(drop=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(df[["sample", "AEI_percent", "AG_mismatches", "A_coverage",
              "signal_to_noise"]].to_string(index=False))


if __name__ == "__main__":
    main()
