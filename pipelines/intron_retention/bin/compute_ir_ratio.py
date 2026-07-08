#!/usr/bin/env python3
"""
compute_ir_ratio.py -- combine intronic + exonic featureCounts into per-intron IR ratios.

IR ratio (coverage-normalized, IRFinder-like):
    intron_density = intron_read_count / intron_length      (reads per bp of pure intron)
    exon_density   = host_gene_exon_count / host_gene_exon_length
    IR_ratio       = intron_density / (intron_density + exon_density)

IR_ratio in [0,1]: 0 = intron fully spliced out; ->1 = intron retained at the
same or higher depth than its host gene's exons. This is the per-intron,
per-sample value the modeling layer needs.

Inputs (featureCounts default text output, one column per sample BAM):
  --intron-counts  featureCounts on introns.saf  (meta-feature = intron_id)
  --exon-counts    featureCounts on exons.saf     (meta-feature = gene_id)
  --map            intron2gene.tsv from make_intron_saf.py
  --run-accession  sample key (this run's ENA run accession) -- when a single
                   sample is processed per invocation. If featureCounts holds
                   multiple BAMs, columns are auto-derived from BAM basenames.
  --cohort         cohort label for the contract's 2nd column

Outputs:
  --out-long     per-intron long table (intron_id, gene_id, IR_ratio, counts, lengths)
  --out-wide     tidy-wide parquet: rows=run_accession, cols=intron_id, values=IR_ratio
                 (first two cols run_accession, cohort) -- the hand-off contract format
  --out-summary  per-sample summary (n_introns_evaluated, median/mean IR, n_high_IR)
"""
import argparse
import os
import sys
import pandas as pd


def read_featurecounts(path):
    """Read a featureCounts output. Returns (df_counts, length_series, sample_cols).
    df_counts indexed by Geneid; sample columns are the trailing BAM columns."""
    df = pd.read_csv(path, sep="\t", comment="#")
    meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    df = df.set_index("Geneid")
    return df, df["Length"], sample_cols


def clean_sample_name(col, override):
    if override:
        return override
    base = os.path.basename(col)
    for suf in (".bam", ".sorted", ".markdup", ".Aligned", ".genome",
                ".Aligned.sortedByCoord.out", ".bai"):
        base = base.replace(suf, "")
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--intron-counts", required=True)
    ap.add_argument("--exon-counts", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--run-accession", default=None,
                    help="sample key override (single-sample invocation)")
    ap.add_argument("--cohort", default="NA")
    ap.add_argument("--min-gene-exon-count", type=int, default=20,
                    help="require host gene to have >= this many exonic reads to trust IR [20]")
    ap.add_argument("--high-ir-threshold", type=float, default=0.1,
                    help="IR ratio above which an intron is called 'retained' for the summary [0.1]")
    ap.add_argument("--out-long", required=True)
    ap.add_argument("--out-wide", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    imap = pd.read_csv(args.map, sep="\t")  # intron_id, gene_id, intron_length
    imap = imap.set_index("intron_id")

    idf, ilen, isamps = read_featurecounts(args.intron_counts)
    edf, elen, esamps = read_featurecounts(args.exon_counts)

    # Both featureCounts tables are run with `-f` (feature-level), so a single
    # logical meta-feature can span MANY rows:
    #   * exon table  — one row per exon, so a gene appears on ~6 rows on average
    #     (360k exon rows for 63k genes). IR needs per-GENE exonic signal.
    #   * intron table — make_intron_saf.py fragments each pure-intronic region
    #     into the disjoint sub-intervals left after masking overlapping exons, so
    #     one intron_id (e.g. GENE__intron_5) appears on several rows.
    # Collapse BOTH to one row per meta-feature id: SUM the read counts and SUM the
    # feature lengths. Without this, the `imap.join(intron_count)` and
    # `d.join(gene_count, on="gene_id")` below join each key against every duplicate
    # row — a cartesian blow-up that OOM-kills the process (and, more subtly,
    # double-counts) on a real genome-wide sample.
    def collapse_by_id(df, length_series):
        if not df.index.has_duplicates:
            return df, length_series
        agg = df.groupby(level=0).sum(numeric_only=True)
        return agg, agg["Length"]

    idf, ilen = collapse_by_id(idf, ilen)
    edf, elen = collapse_by_id(edf, elen)

    # Pair intron and exon columns by cleaned sample NAME (not by position), so
    # the two featureCounts runs align correctly even if column order differs.
    single = (len(isamps) == 1 and len(esamps) == 1)
    ekey = {clean_sample_name(c, None): c for c in esamps}
    pairs = []  # list of (sample_name, intron_col, exon_col)
    if single:
        # single-sample invocation: honour --run-accession override, pair the lone cols
        samp = clean_sample_name(isamps[0], args.run_accession)
        pairs.append((samp, isamps[0], esamps[0]))
    else:
        unmatched = []
        for icol in isamps:
            sname = clean_sample_name(icol, None)
            ecol = ekey.get(sname)
            if ecol is None:
                unmatched.append(sname)
                continue
            pairs.append((sname, icol, ecol))
        if unmatched:
            sys.stderr.write(
                "[compute_ir_ratio] WARNING: no exon column matched intron sample(s) "
                f"{unmatched}; these are skipped\n")
        if not pairs:
            sys.stderr.write("[compute_ir_ratio] ERROR: no intron/exon sample columns "
                             "matched by name\n")
            sys.exit(1)

    long_frames = []
    wide_rows = {}
    summary_rows = []

    for samp, icol, ecol in pairs:

        intron_count = idf[icol]
        gene_count = edf[ecol]
        gene_exlen = elen  # summed exonic length per gene

        # align introns to host genes
        d = imap.join(intron_count.rename("intron_count"), how="inner")
        d.index.name = "intron_id"  # join drops it when index names differ
        d["intron_length"] = d["intron_length"].astype(float)
        d = d.join(gene_count.rename("gene_exon_count"), on="gene_id")
        d = d.join(gene_exlen.rename("gene_exon_length"), on="gene_id")
        d = d.dropna(subset=["gene_exon_count", "gene_exon_length"])

        intron_density = d["intron_count"] / d["intron_length"]
        exon_density = d["gene_exon_count"] / d["gene_exon_length"]
        denom = intron_density + exon_density
        ir = intron_density / denom
        ir[denom == 0] = pd.NA  # gene silent -> undefined
        # require enough host-gene exonic signal to trust the ratio
        ir[d["gene_exon_count"] < args.min_gene_exon_count] = pd.NA

        d = d.assign(run_accession=samp, cohort=args.cohort,
                     intron_density=intron_density, exon_density=exon_density,
                     IR_ratio=ir)
        long_frames.append(d.reset_index()[[
            "run_accession", "cohort", "intron_id", "gene_id",
            "intron_count", "intron_length", "gene_exon_count",
            "gene_exon_length", "IR_ratio"]])

        wide_rows[samp] = ir

        valid = ir.dropna()
        summary_rows.append({
            "run_accession": samp,
            "cohort": args.cohort,
            "n_introns_total": int(len(ir)),
            "n_introns_evaluated": int(valid.shape[0]),
            "median_IR": float(valid.median()) if len(valid) else float("nan"),
            "mean_IR": float(valid.mean()) if len(valid) else float("nan"),
            f"n_IR_gt_{args.high_ir_threshold}": int((valid > args.high_ir_threshold).sum()),
        })

    long_df = pd.concat(long_frames, ignore_index=True)
    long_df.to_csv(args.out_long, sep="\t", index=False)

    # tidy-wide: rows = samples, cols = intron_id, values = IR_ratio
    wide = pd.DataFrame(wide_rows).T  # samples x introns
    wide.index.name = "run_accession"
    wide = wide.reset_index()
    wide.insert(1, "cohort", args.cohort)
    # numeric feature columns as float32 to keep the matrix small
    feat_cols = [c for c in wide.columns if c not in ("run_accession", "cohort")]
    # cast column-by-column: a single assignment `wide[feat_cols] = ...` fails with
    # "Columns must be same length as key" if any intron_id is duplicated across the
    # feature columns (defensive — the collapse above should already make them unique).
    for c in feat_cols:
        wide[c] = wide[c].astype("float32")
    try:
        wide.to_parquet(args.out_wide, index=False, compression="gzip")
    except Exception as e:  # pragma: no cover - pyarrow missing
        sys.stderr.write(f"[compute_ir_ratio] parquet failed ({e}); writing CSV\n")
        wide.to_csv(args.out_wide.replace(".parquet", ".csv"), index=False)

    pd.DataFrame(summary_rows).to_csv(args.out_summary, sep="\t", index=False)
    sys.stderr.write(
        f"[compute_ir_ratio] {len(summary_rows)} sample(s); "
        f"{long_df['intron_id'].nunique()} introns; wide matrix {wide.shape}\n"
    )


if __name__ == "__main__":
    main()
