#!/usr/bin/env python3
"""
ir_antigen_features.py — CLI wrapper that turns the pipeline's intron-retention
output into the two NAMED per-sample features (retained_intron_load,
ir_neoantigen_burden) via analysis/differentiated/intron_retention.py.

Consumes:
  --ir-matrix    intron_retention.parquet (contract tidy-wide) OR --ir-long *.ir_long.tsv
  --intron-saf   introns.saf              (from make_intron_saf.py, coordinates)
  --genome       GRCh38 primary_assembly FASTA (indexed .fai for the pilot)
  --hla-table    hla_typing.parquet       (per-sample HLA-I alleles)
Emits:
  --out          feature parquet keyed on (run_accession, cohort) with:
                 retained_intron_load, n_introns_evaluated, retained_intron_fraction,
                 retained_intron_load_weighted, retained_intron_load_cohortz,
                 ir_neoantigen_burden, n_candidate_peptides, n_retained_introns_used

This is the process body the intron_retention.nf subworkflow calls; the feature
logic itself lives (and is unit-tested) in intron_retention.py.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

# make the module importable whether run from repo root or from bin/
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))            # analysis/differentiated
sys.path.insert(0, str(_HERE.parents[2] / "antigen_core"))  # analysis/antigen_core
import intron_retention as ir  # noqa: E402


def _read_ir(args) -> pd.DataFrame:
    if args.ir_matrix:
        return pd.read_parquet(args.ir_matrix)
    frames = [pd.read_csv(f, sep="\t") for f in args.ir_long]
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ir-matrix", help="intron_retention.parquet (tidy-wide)")
    src.add_argument("--ir-long", nargs="+", help="per-sample *.ir_long.tsv files")
    ap.add_argument("--intron-saf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--hla-table", required=True)
    ap.add_argument("--threshold", type=float, default=ir.IR_RETAINED_THRESHOLD)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ir_data = _read_ir(args)
    hla_table = pd.read_parquet(args.hla_table)
    feat = ir.build_ir_features(
        ir_data, args.intron_saf, args.genome, hla_table, threshold=args.threshold)

    out = Path(args.out)
    try:
        feat.to_parquet(out, index=False)
    except Exception as e:  # pragma: no cover
        sys.stderr.write(f"[ir_antigen_features] parquet failed ({e}); CSV fallback\n")
        feat.to_csv(str(out).replace(".parquet", ".csv"), index=False)
    sys.stderr.write(
        f"[ir_antigen_features] {feat.shape[0]} samples; "
        f"retained_intron_load median="
        f"{feat['retained_intron_load'].median()}, "
        f"ir_neoantigen_burden non-NA="
        f"{feat['ir_neoantigen_burden'].notna().sum()}\n")


if __name__ == "__main__":
    main()
