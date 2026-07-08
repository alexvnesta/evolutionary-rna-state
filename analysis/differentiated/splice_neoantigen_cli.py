#!/usr/bin/env python
"""
splice_neoantigen_cli.py — thin CLI shim for splicing_neoantigen.nf.

Single-sources the splicing-neoantigen ALGORITHM with the unit-tested
``analysis/differentiated/splicing_neoantigen.py`` module: each Nextflow
process calls one subcommand here rather than duplicating logic. Subcommands
mirror the subworkflow steps:

    call-neojunctions   SJ.out.tab            -> neojunctions.tsv  (SNAF gate)
    translate           neojunctions + genome -> peptides.txt      (SNAF 3-frame)
    score-burden        peptides + HLA table  -> burden.tsv        (shared engine)
    merge-burden        per-sample burden tsv -> feature parquet

Kept deliberately small; all science lives in the importable module.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# repo root on path so `analysis...` imports resolve when run from Nextflow work dir
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd

from analysis.differentiated.splicing_neoantigen import (
    Junction, FastaSeq, ExonIndex,
    parse_star_sj, call_neojunctions,
    peptides_from_neojunctions, splice_neoantigen_burden,
    junctions_from_frame, KS, DEFAULT_FLANK, T_MIN, N_MAX,
)

__version__ = "1.0.0-snaf-algorithm-port"


def _load_hla_alleles(hla_table: str, run_accession: str) -> list[str]:
    """Pull one sample's 6 HLA-I alleles from the hla_typing table."""
    df = pd.read_parquet(hla_table) if str(hla_table).endswith(".parquet") \
        else pd.read_csv(hla_table, sep=None, engine="python")
    row = df[df["run_accession"] == run_accession]
    if row.empty:
        return []
    cols = ["HLA_A_1", "HLA_A_2", "HLA_B_1", "HLA_B_2", "HLA_C_1", "HLA_C_2"]
    return [str(row.iloc[0][c]) for c in cols if c in df.columns and pd.notna(row.iloc[0][c])]


def cmd_call_neojunctions(a):
    junctions = parse_star_sj(a.sj, min_reads=1)
    if a.normal_ref:
        # normal_ref: tidy tsv chrom,start,end,strand,mean_count -> subtract means
        nr = pd.read_csv(a.normal_ref, sep="\t")
        key = {(r.chrom, int(r.start), int(r.end), r.strand): float(r.mean_count)
               for r in nr.itertuples()}
        for j in junctions:
            j.normal_mean = key.get((j.chrom, j.start, j.end, j.strand), 0.0)
    neo = call_neojunctions(junctions, t_min=a.t_min, n_max=a.n_max)
    with open(a.out, "w") as fh:
        fh.write("chrom\tstart\tend\tstrand\tcount\tnormal_mean\n")
        for j in neo:
            fh.write(f"{j.chrom}\t{j.start}\t{j.end}\t{j.strand}\t{j.count}\t{j.normal_mean}\n")


def cmd_translate(a):
    df = pd.read_csv(a.neojunctions, sep="\t")
    junctions = junctions_from_frame(
        df, count_col="count",
        normal_mean_col="normal_mean" if "normal_mean" in df.columns else None)
    chroms = {j.chrom for j in junctions}
    fasta = FastaSeq(a.fasta)
    exon_index = ExonIndex.from_gtf(a.gtf, chroms=chroms) if a.gtf else None
    ks = tuple(int(x) for x in a.ks.split(","))
    peps = peptides_from_neojunctions(junctions, fasta, exon_index=exon_index,
                                      ks=ks, flank=a.flank)
    with open(a.out, "w") as fh:
        fh.write("\n".join(peps) + ("\n" if peps else ""))


def cmd_score_burden(a):
    peps = [l.strip() for l in open(a.peptides) if l.strip()]
    alleles = _load_hla_alleles(a.hla_table, a.run_accession)
    # peptides are already derived; score directly through the shared engine
    from analysis.antigen_core.mhc_binding import count_binders
    burden = count_binders(peps, alleles, rank_threshold=a.rank)
    with open(a.out, "w") as fh:
        fh.write("run_accession\tcohort\tsplice_neoantigen_burden\n")
        fh.write(f"{a.run_accession}\t{a.cohort}\t{int(burden)}\n")


def cmd_merge_burden(a):
    frames = [pd.read_csv(p, sep="\t") for p in a.burdens]
    out = pd.concat(frames, ignore_index=True)
    out = out[["run_accession", "cohort", "splice_neoantigen_burden"]]
    out["splice_neoantigen_burden"] = out["splice_neoantigen_burden"].astype(int)
    out.to_parquet(a.out, index=False)


def main():
    p = argparse.ArgumentParser(prog="splice_neoantigen_cli.py")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("call-neojunctions")
    s1.add_argument("--sj", required=True)
    s1.add_argument("--t-min", dest="t_min", type=int, default=T_MIN)
    s1.add_argument("--n-max", dest="n_max", type=int, default=N_MAX)
    s1.add_argument("--normal-ref", dest="normal_ref", default=None)
    s1.add_argument("--out", required=True)
    s1.set_defaults(func=cmd_call_neojunctions)

    s2 = sub.add_parser("translate")
    s2.add_argument("--neojunctions", required=True)
    s2.add_argument("--fasta", required=True)
    s2.add_argument("--gtf", default=None)
    s2.add_argument("--ks", default=",".join(str(k) for k in KS))
    s2.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    s2.add_argument("--out", required=True)
    s2.set_defaults(func=cmd_translate)

    s3 = sub.add_parser("score-burden")
    s3.add_argument("--peptides", required=True)
    s3.add_argument("--hla-table", dest="hla_table", required=True)
    s3.add_argument("--run-accession", dest="run_accession", required=True)
    s3.add_argument("--cohort", required=True)
    s3.add_argument("--rank", type=float, default=2.0)
    s3.add_argument("--out", required=True)
    s3.set_defaults(func=cmd_score_burden)

    s4 = sub.add_parser("merge-burden")
    s4.add_argument("--burdens", nargs="+", required=True)
    s4.add_argument("--out", required=True)
    s4.set_defaults(func=cmd_merge_burden)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
