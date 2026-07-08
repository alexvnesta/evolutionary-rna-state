#!/usr/bin/env python3
"""
make_alu_bed.py -- build an Alu-only BED6 from the UCSC hg38 RepeatMasker table.

Input : rmsk.txt(.gz) from
        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
        (columns include: swScore genoName genoStart genoEnd strand repName
         repClass repFamily -- repClass 'SINE', repFamily 'Alu')
Output: BED6 (chrom, start, end, repName, swScore, strand), Alu family only.

Contig naming: UCSC rmsk uses 'chr1..chrX,chrY,chrM'. GENCODE primary_assembly
FASTA also uses 'chr'-prefixed names, so no renaming is required for this repo's
reference. A --keep-standard flag drops _alt/_random/_fix scaffolds.

Usage:
    make_alu_bed.py --rmsk rmsk.txt.gz --out alu.hg38.bed6 [--keep-standard]
"""
import argparse
import gzip
import sys

STD = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep-standard", action="store_true",
                    help="keep only chr1..22,X,Y,M (drop _alt/_random/_fix)")
    args = ap.parse_args()

    n_in = n_out = 0
    with opener(args.rmsk) as fh, open(args.out, "w") as o:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            # rmsk.txt has a leading 'bin' column, so fields shift by +1 vs rmsk.sql
            # Layout: bin swScore milliDiv milliDel milliIns genoName genoStart
            #         genoEnd genoLeft strand repName repClass repFamily ...
            if len(f) < 13:
                continue
            genoName, genoStart, genoEnd = f[5], f[6], f[7]
            strand, repName, repClass, repFamily = f[9], f[10], f[11], f[12]
            n_in += 1
            if repClass != "SINE" or repFamily != "Alu":
                continue
            if args.keep_standard and genoName not in STD:
                continue
            strand = "+" if strand == "+" else "-"
            swScore = f[1]
            o.write(f"{genoName}\t{genoStart}\t{genoEnd}\t{repName}\t{swScore}\t{strand}\n")
            n_out += 1
    sys.stderr.write(f"[make_alu_bed] {n_out} Alu intervals written from {n_in} rmsk rows\n")


if __name__ == "__main__":
    main()
