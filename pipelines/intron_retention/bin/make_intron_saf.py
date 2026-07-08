#!/usr/bin/env python3
"""
make_intron_saf.py  --  derive pure-intronic and exonic intervals from a GENCODE GTF
                        and emit two featureCounts SAF files.

Intron-retention (IR) quantification, featureCounts-based (arm64-native) approach.

Definitions (IRFinder-like, but computed from annotation only, no source build):
  * For every gene, take the UNION of all exons across all transcripts (the
    "exonic" footprint). The gaps between consecutive union-exon blocks are the
    candidate intron blocks of that gene.
  * Each candidate intron block is then MASKED against the genome-wide union of
    ALL exons (every gene, both strands). Any sub-region of a candidate intron
    that overlaps some other feature's exon is removed. What survives is
    "pure intronic" sequence -- present in the pre-mRNA, never exonic in any
    annotated isoform of any gene. Reads there are retained-intron signal, not
    spliced-in exon signal from an overlapping/antisense gene.

Outputs (SAF = GeneID, Chr, Start, End, Strand; 1-based inclusive, featureCounts):
  introns.saf : one row per pure-intronic sub-interval. Rows are grouped into
                meta-features by GeneID = "<gene_id>__intron_<n>" so featureCounts
                sums a fragmented intron back into a single count + summed Length.
  exons.saf   : one row per union-exon block, GeneID = gene_id, so featureCounts
                yields per-gene total exonic count + total exonic Length.
  intron2gene.tsv : intron_id <tab> gene_id <tab> intron_length  (host-gene map)

The two SAF files feed two featureCounts runs; compute_ir_ratio.py combines them.
Pure stdlib -- no pandas/pyranges -- so the conda env stays minimal on arm64.
"""
import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_gene_id(attr):
    # attribute field: gene_id "ENSG..."; transcript_id "..."; ...
    # fast targeted parse (avoid splitting the whole field)
    key = 'gene_id "'
    i = attr.find(key)
    if i < 0:
        return None
    i += len(key)
    j = attr.find('"', i)
    return attr[i:j]


def merge_intervals(ivs):
    """Merge a list of (start, end) 1-based inclusive intervals. Returns sorted merged list."""
    if not ivs:
        return []
    ivs.sort()
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:  # touching or overlapping -> merge
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def subtract_mask(block_s, block_e, mask_sorted):
    """Subtract merged mask intervals (sorted, non-overlapping) from [block_s, block_e].
    Returns list of surviving (start, end) 1-based inclusive sub-intervals."""
    survivors = []
    cur = block_s
    for ms, me in mask_sorted:
        if me < cur:
            continue
        if ms > block_e:
            break
        # mask overlaps [cur, block_e]
        if ms > cur:
            survivors.append((cur, min(ms - 1, block_e)))
        cur = max(cur, me + 1)
        if cur > block_e:
            break
    if cur <= block_e:
        survivors.append((cur, block_e))
    return survivors


def overlapping_slice(mask_sorted, s, e):
    """Return the sub-list of mask intervals that could overlap [s, e].
    mask_sorted is sorted by start. Linear scan bounded (chromosome-local)."""
    # simple: filter (chromosome-scale lists are fine for a single-machine run)
    return [(ms, me) for ms, me in mask_sorted if not (me < s or ms > e)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="GENCODE/Ensembl GTF (.gtf or .gtf.gz)")
    ap.add_argument("--out-introns", required=True, help="output introns.saf")
    ap.add_argument("--out-exons", required=True, help="output exons.saf")
    ap.add_argument("--out-map", required=True, help="output intron2gene.tsv")
    ap.add_argument("--min-intron-len", type=int, default=50,
                    help="drop pure-intronic sub-intervals shorter than this (bp) [50]")
    args = ap.parse_args()

    # gene_id -> chrom, strand, list of exon (start,end)
    gene_exons = defaultdict(list)
    gene_meta = {}
    # chrom -> list of exon (start,end)   (genome-wide exon mask, strand-agnostic)
    chrom_exons = defaultdict(list)

    n_exon = 0
    with open_maybe_gz(args.gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom = f[0]
            start = int(f[3])
            end = int(f[4])
            strand = f[6]
            gid = parse_gene_id(f[8])
            if gid is None:
                continue
            gene_exons[gid].append((start, end))
            chrom_exons[chrom].append((start, end))
            if gid not in gene_meta:
                gene_meta[gid] = (chrom, strand)
            n_exon += 1
    sys.stderr.write(f"[make_intron_saf] parsed {n_exon} exon lines, {len(gene_exons)} genes\n")

    # genome-wide merged exon mask per chromosome
    chrom_mask = {c: merge_intervals(ivs) for c, ivs in chrom_exons.items()}
    sys.stderr.write("[make_intron_saf] built genome-wide exon mask\n")

    fi = open(args.out_introns, "w")
    fe = open(args.out_exons, "w")
    fm = open(args.out_map, "w")
    fi.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    fe.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    fm.write("intron_id\tgene_id\tintron_length\n")

    n_introns = 0
    n_intron_rows = 0
    n_genes_with_intron = 0
    for gid, exons in gene_exons.items():
        chrom, strand = gene_meta[gid]
        union = merge_intervals(exons)  # union-exon blocks, sorted
        # write exon rows (per-gene meta-feature)
        for s, e in union:
            fe.write(f"{gid}\t{chrom}\t{s}\t{e}\t{strand}\n")
        if len(union) < 2:
            continue  # single-exon gene -> no introns
        mask = chrom_mask[chrom]
        # candidate intron blocks = gaps between consecutive union-exon blocks
        gaps = []
        for k in range(len(union) - 1):
            gs = union[k][1] + 1
            ge = union[k + 1][0] - 1
            if ge >= gs:
                gaps.append((gs, ge))
        if not gaps:
            continue
        # strand-aware intron numbering: 5'->3'
        ordered = gaps if strand != "-" else list(reversed(gaps))
        gene_had_intron = False
        for idx, (gs, ge) in enumerate(ordered, start=1):
            # mask against genome-wide exons (removes overlapping/antisense exons)
            local_mask = overlapping_slice(mask, gs, ge)
            survivors = subtract_mask(gs, ge, local_mask)
            survivors = [(s, e) for s, e in survivors if (e - s + 1) >= args.min_intron_len]
            if not survivors:
                continue
            intron_id = f"{gid}__intron_{idx}"
            ilen = sum(e - s + 1 for s, e in survivors)
            for s, e in survivors:
                fi.write(f"{intron_id}\t{chrom}\t{s}\t{e}\t{strand}\n")
                n_intron_rows += 1
            fm.write(f"{intron_id}\t{gid}\t{ilen}\n")
            n_introns += 1
            gene_had_intron = True
        if gene_had_intron:
            n_genes_with_intron += 1

    fi.close()
    fe.close()
    fm.close()
    sys.stderr.write(
        f"[make_intron_saf] wrote {n_introns} pure-intronic meta-features "
        f"({n_intron_rows} rows) across {n_genes_with_intron} genes\n"
    )


if __name__ == "__main__":
    main()
