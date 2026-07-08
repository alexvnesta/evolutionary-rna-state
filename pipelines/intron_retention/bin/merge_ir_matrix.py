#!/usr/bin/env python3
"""
merge_ir_matrix.py -- assemble per-sample IR long tables into one cohort matrix.

Output matches docs/HANDOFF_CONTRACT.md:
  intron_retention.parquet -- rows = run_accession, first two cols
                              (run_accession, cohort), remaining cols = intron_id,
                              values = IR ratio; NA where a sample lacked the intron.
Also concatenates the per-sample summaries into one table.
"""
import argparse
import sys
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--long", nargs="+", required=True, help="per-sample *.ir_long.tsv")
    ap.add_argument("--summaries", nargs="+", required=True, help="per-sample *.ir_summary.tsv")
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    frames = []
    for f in args.long:
        d = pd.read_csv(f, sep="\t", usecols=["run_accession", "cohort", "intron_id", "IR_ratio"])
        frames.append(d)
    allrows = pd.concat(frames, ignore_index=True)

    # pivot to samples x introns
    wide = allrows.pivot_table(index=["run_accession", "cohort"],
                               columns="intron_id", values="IR_ratio",
                               aggfunc="first")
    wide = wide.reset_index()
    # Down-cast the intron block to float32 in ONE block op (not element/column-wise).
    # With ~290k intron columns, `wide[feat_cols] = wide[feat_cols].astype(...)` is
    # O(ncols) column re-inserts and stalls for minutes; drop+astype+concat is O(1)
    # in column count. feat_cols are unique (pivot dedupes intron_id).
    id_cols = ["run_accession", "cohort"]
    feat_cols = [c for c in wide.columns if c not in id_cols]
    feat = wide.drop(columns=id_cols).astype("float32")
    wide = pd.concat([wide[id_cols], feat], axis=1)
    try:
        wide.to_parquet(args.out_matrix, index=False, compression="gzip")
    except Exception as e:  # pragma: no cover
        sys.stderr.write(f"[merge_ir_matrix] parquet failed ({e}); CSV fallback\n")
        wide.to_csv(args.out_matrix.replace(".parquet", ".csv"), index=False)

    sumframes = [pd.read_csv(f, sep="\t") for f in args.summaries]
    pd.concat(sumframes, ignore_index=True).to_csv(args.out_summary, sep="\t", index=False)
    sys.stderr.write(f"[merge_ir_matrix] cohort matrix {wide.shape} "
                     f"({len(feat_cols)} introns, {wide.shape[0]} samples)\n")


if __name__ == "__main__":
    main()
