#!/usr/bin/env python3
"""
compute_aei.py -- Alu Editing Index (AEI) for a single RNA-seq BAM.

The AEI (Roth, Levanon et al., Nat Methods 2019) is a robust, cohort-comparable
summary of global A-to-I editing: the number of A>G mismatches divided by the
total coverage of adenosines, aggregated over *all* Alu elements genome-wide.
Because it pools millions of Alu positions it is stable at modest depth and is
far less sensitive to per-site noise than counting individual edited sites.

Definition (strand-aware):
    For each Alu interval, orient reads to the sense strand of the Alu.
    Sense A positions:  edited signal = #G  , reference signal = #A
    (on '-' strand Alus we look at T/C on the genomic forward strand, which is
     A/G on the Alu sense strand -- handled by strand flipping below).

    AEI (%) = 100 * sum(mismatch_A>G) / sum(A>G-informative coverage over A refs)

This script uses a genome pileup restricted to Alu intervals via pysam and the
reference FASTA to know which genomic base is a reference 'A' (sense) / 'T'
(antisense). It reports the index plus the raw numerator/denominator and the
per-mismatch-type breakdown (a QC control: A>G should dominate; other
mismatch types estimate the noise floor).

Usage:
    compute_aei.py --bam s.bam --fasta genome.fa --alu alu.bed \
        --sample S1 --out S1.aei.tsv [--min-baseq 25 --min-mapq 255 \
        --snp-bed dbsnp.bed.gz]

Notes for arm64: pysam wheels/conda builds are native osx-arm64; no Docker.
"""
import argparse
import sys
import pysam

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
# On the Alu SENSE strand we score A->G editing. If the Alu is on the genomic
# '-' strand, the sense 'A' corresponds to a genomic 'T', and sense 'G' to a
# genomic 'C'. We therefore flip observed bases for '-' Alus.


def load_intervals(bed):
    ivs = []
    with open(bed) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            strand = f[5] if len(f) > 5 else "+"
            ivs.append((chrom, start, end, strand))
    return ivs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--alu", required=True, help="Alu intervals BED (0-based, +/- strand col 6)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-baseq", type=int, default=25)
    ap.add_argument("--min-mapq", type=int, default=255,
                    help="STAR uniquely-mapped reads have MAPQ 255; default keeps only those")
    ap.add_argument("--snp-bed", default=None,
                    help="Optional bgzipped+tabixed BED of known SNPs to mask (dbSNP)")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)
    snp = None
    if args.snp_bed:
        snp = pysam.TabixFile(args.snp_bed)
        snp_contigs = set(snp.contigs)

    # 12 directed mismatch types on the Alu sense strand.
    types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
             "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
    counts = {t: 0 for t in types}
    ref_cov = {b: 0 for b in "ACGT"}   # coverage over reference bases (sense)

    bam_contigs = set(bam.references)
    fasta_contigs = set(fasta.references)

    for chrom, start, end, strand in load_intervals(args.alu):
        if chrom not in bam_contigs or chrom not in fasta_contigs:
            continue
        refseq = fasta.fetch(chrom, start, end).upper()
        for col in bam.pileup(chrom, start, end, truncate=True,
                              stepper="samtools", min_base_quality=args.min_baseq,
                              ignore_orphans=False):
            pos = col.reference_pos
            rb = refseq[pos - start]
            if rb not in "ACGT":
                continue
            # dbSNP masking: skip positions overlapping a known SNP
            if snp is not None and chrom in snp_contigs:
                if any(True for _ in snp.fetch(chrom, pos, pos + 1)):
                    continue
            # Sense-strand reference base and observed-base flip for '-' Alus
            sense_ref = rb if strand == "+" else COMPLEMENT[rb]
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality < args.min_mapq:
                    continue
                ob = aln.query_sequence[pr.query_position].upper()
                if ob not in "ACGT":
                    continue
                sense_ob = ob if strand == "+" else COMPLEMENT[ob]
                ref_cov[sense_ref] += 1
                if sense_ob != sense_ref:
                    counts[f"{sense_ref}>{sense_ob}"] += 1

    # AEI = A>G editing rate over adenosine coverage (sense strand), as percent.
    a_cov = ref_cov["A"]
    ag = counts["A>G"]
    aei = 100.0 * ag / a_cov if a_cov else 0.0
    # Noise floor proxy: mean of the non-A>G mismatch rates over their ref coverage
    def rate(mm):
        r = mm[0]
        return counts[mm] / ref_cov[r] if ref_cov[r] else 0.0
    control_types = [t for t in types if t != "A>G"]
    noise = 100.0 * sum(rate(t) for t in control_types) / len(control_types)

    with open(args.out, "w") as o:
        o.write("sample\tAEI_percent\tAG_mismatches\tA_coverage\t"
                "signal_to_noise\tnoise_floor_percent\t"
                + "\t".join(f"cov_{b}" for b in "ACGT") + "\t"
                + "\t".join(f"n_{t.replace('>','to')}" for t in types) + "\n")
        s2n = (aei / noise) if noise else float("inf")
        o.write(f"{args.sample}\t{aei:.6f}\t{ag}\t{a_cov}\t{s2n:.3f}\t{noise:.6f}\t"
                + "\t".join(str(ref_cov[b]) for b in "ACGT") + "\t"
                + "\t".join(str(counts[t]) for t in types) + "\n")
    sys.stderr.write(f"[AEI] {args.sample}: AEI={aei:.4f}%  A>G={ag}  A_cov={a_cov}  S/N={s2n:.2f}\n")


if __name__ == "__main__":
    main()
