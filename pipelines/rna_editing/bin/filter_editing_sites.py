#!/usr/bin/env python3
"""
filter_editing_sites.py -- extract A-to-I sites from a JACUSA2 call-1 table.

JACUSA2 `call-1` writes a header line beginning with '#contig' and then one row
per candidate position. The per-sample base-count column ("bases11") is a
comma-separated A,C,G,T count. A-to-I editing appears as:
    sense '+' : reference A, alternate G   -> freq = G / (A+G)
    sense '-' : reference T, alternate C   -> freq = C / (C+T)  (genomic strand)

We keep sites where coverage, edited-read count and editing frequency pass the
given thresholds, and (optionally) drop positions listed in a bgzipped+tabixed
dbSNP BED. Output is a tidy TSV of edited sites with frequency and depth.

This parser is defensive about column order: it locates 'bases11' and 'strand'
from the header, and falls back to fixed BED6+counts positions if absent.
"""
import argparse
import sys

try:
    import pysam
except ImportError:
    pysam = None


def parse_counts(field):
    # "A,C,G,T" integer counts
    a, c, g, t = (int(x) for x in field.split(","))
    return a, c, g, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jacusa", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--min-edit-freq", type=float, default=0.10)
    ap.add_argument("--min-edit-reads", type=int, default=3)
    ap.add_argument("--snp-bed", default=None)
    args = ap.parse_args()

    snp = None
    if args.snp_bed and pysam is not None:
        snp = pysam.TabixFile(args.snp_bed)
        snp_contigs = set(snp.contigs)

    hdr_cols = None
    bidx = sidx = ridx = None
    n_kept = 0
    with open(args.jacusa) as fh, open(args.out, "w") as o:
        o.write("chrom\tpos\tstrand\tref\talt\tedit_freq\tedited_reads\tcoverage\n")
        for ln in fh:
            if ln.startswith("#"):
                if ln.startswith("#contig") or ln.lower().startswith("#chrom"):
                    hdr_cols = ln.lstrip("#").rstrip("\n").split("\t")
                    for i, c in enumerate(hdr_cols):
                        if c.startswith("bases"):
                            bidx = i
                        if c == "strand":
                            sidx = i
                        if c == "ref":
                            ridx = i
                continue
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            chrom = f[0]
            end = int(f[2])          # BED end == 1-based position
            strand = f[sidx] if sidx is not None and sidx < len(f) else f[5]
            counts_field = f[bidx] if bidx is not None and bidx < len(f) else f[6]
            refbase = (f[ridx].upper() if ridx is not None and ridx < len(f) else "")
            try:
                a, c, g, t = parse_counts(counts_field)
            except (ValueError, IndexError):
                continue
            cov = a + c + g + t
            if cov < args.min_cov:
                continue
            # dbSNP masking
            if snp is not None and chrom in snp_contigs:
                if any(True for _ in snp.fetch(chrom, end - 1, end)):
                    continue
            # Decide A-to-I orientation.
            #   stranded '+' : ref A, edit -> G ;  '-' : ref T, edit -> C
            #   unstranded '.' : JACUSA2 reports the genomic reference base, so use
            #     it directly — ref A => A>G, ref T => T>C. This is the only way to
            #     score A-to-I when the library strandedness is unknown/unset;
            #     without it every '.' site is dropped and 0 edits are ever kept.
            if strand == "+":
                ref, alt, edited, denom = "A", "G", g, a + g
            elif strand == "-":
                ref, alt, edited, denom = "T", "C", c, c + t
            elif refbase == "A":
                ref, alt, edited, denom = "A", "G", g, a + g
            elif refbase == "T":
                ref, alt, edited, denom = "T", "C", c, c + t
            else:
                continue
            if denom == 0:
                continue
            freq = edited / denom
            if edited >= args.min_edit_reads and freq >= args.min_edit_freq:
                out_strand = strand if strand in ("+", "-") else ("+" if ref == "A" else "-")
                o.write(f"{chrom}\t{end}\t{out_strand}\t{ref}\t{alt}\t{freq:.4f}\t{edited}\t{cov}\n")
                n_kept += 1
    sys.stderr.write(f"[filter_editing_sites] {n_kept} A-to-I sites kept\n")


if __name__ == "__main__":
    main()
