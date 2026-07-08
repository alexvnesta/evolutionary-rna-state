#!/usr/bin/env python
"""
bin/build_te_antigen_table.py — per-sample TE/ERV antigen burden -> cohort
tidy feature matrix (contract format). Called by the TE_ANTIGEN Nextflow
subworkflow's BUILD_TE_ANTIGEN_TABLE process, one-shot per cohort.

Inputs (all produced upstream by the pipeline session, per HANDOFF_CONTRACT):
  --te-locus      te_locus.parquet   Telescope per-locus counts
                                     (rows=run_accession, cols=locus_id).
  --annotation    RepeatMasker/GENCODE-TE table: locus_id, repeat_class,
                  chrom, start, end, strand.
  --genome        genome FASTA (+ .fai) the pipeline aligned against — used to
                  extract expressed-locus sequences (pysam).
  --hla           hla_typing.parquet from analysis/antigen_core/hla_typing.py
                  (run_accession + 6 HLA-I allele columns).
  --out           output te_antigen.parquet.

For each sample: pick expressed loci, extract their sequences ONCE, then call
analysis/differentiated/te_antigen.te_antigen_row through the SHARED engine.
No per-sample values are fabricated — every row comes from real inputs.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# make analysis/ importable (antigen_core + differentiated), whether run from
# the repo or a Nextflow work dir.
_ANALYSIS = Path(__file__).resolve().parents[2]     # .../analysis
for p in (str(_ANALYSIS), str(_ANALYSIS / "differentiated")):
    if p not in sys.path:
        sys.path.insert(0, p)

import pandas as pd  # noqa: E402

import te_antigen as te  # noqa: E402
from antigen_core.hla_typing import ALLELE_COLS  # noqa: E402


def _alleles_for(hla_row: pd.Series) -> list[str]:
    return [str(hla_row[c]) for c in ALLELE_COLS
            if c in hla_row and pd.notna(hla_row[c]) and str(hla_row[c]).strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--te-locus", required=True, help="te_locus.parquet")
    ap.add_argument("--annotation", required=True,
                    help="TE locus annotation (parquet/csv): locus_id, repeat_class, chrom, start, end, strand")
    ap.add_argument("--genome", required=True, help="genome FASTA (+ .fai)")
    ap.add_argument("--hla", required=True, help="hla_typing.parquet")
    ap.add_argument("--min-reads", type=float, default=10.0)
    ap.add_argument("--min-cpm", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    te_locus = pd.read_parquet(args.te_locus)
    ann = (pd.read_parquet(args.annotation)
           if args.annotation.endswith(".parquet")
           else pd.read_csv(args.annotation))
    hla = pd.read_parquet(args.hla).set_index("run_accession")

    locus_cols = [c for c in te_locus.columns if c not in ("run_accession", "cohort")]

    rows = []
    for _, srow in te_locus.iterrows():
        run = srow["run_accession"]
        cohort = srow.get("cohort", "")
        counts = {l: float(srow[l]) for l in locus_cols if pd.notna(srow[l])}
        if run not in hla.index:
            continue                       # no HLA genotype -> skip (NA row)
        alleles = _alleles_for(hla.loc[run])
        active = te.select_expressed_loci(
            counts, min_reads=args.min_reads, min_cpm=args.min_cpm)
        locus_seqs = te.extract_locus_sequences(ann, args.genome, locus_ids=active)
        rows.append(te.te_antigen_row(
            run, cohort, counts, locus_seqs, alleles,
            annotation=ann, min_reads=args.min_reads, min_cpm=args.min_cpm))

    df = te.build_te_antigen_table(rows)
    if args.out.endswith(".parquet"):
        df.to_parquet(args.out, index=False)
    else:
        df.to_csv(args.out, index=False)
    print(f"[build_te_antigen_table] wrote {len(df)} samples -> {args.out}")


if __name__ == "__main__":
    main()
