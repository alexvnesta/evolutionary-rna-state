#!/usr/bin/env python3
"""
compute_aei_fast.py -- Alu Editing Index via a single streaming `samtools mpileup`.

Same definition as the repo's compute_aei.py (Roth/Levanon 2019): AEI = 100 * (A>G
mismatches) / (A coverage), pooled over Alu adenosines, strand-aware. But instead of
1.18M per-interval pysam pileup calls (Python-bound, hours/sample) this streams the
BAM once through samtools mpileup restricted to the Alu BED (-l), then parses the C
output. Strand handled via the Alu BED strand column: '+' Alu scores A(ref)->G; '-'
Alu scores T(ref)->C on the genomic forward strand (= A->G on the Alu sense strand).

Usage: compute_aei_fast.py --bam s.bam --fasta genome.fa --alu alu.bed6 \
         --sample S --min-baseq 25 --min-mapq 60 --out S.aei.tsv
"""
import argparse, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--alu", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-baseq", type=int, default=25)
    ap.add_argument("--min-mapq", type=int, default=60)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # strand lookup per position is too big; instead split the BED by strand and run twice.
    plus_bed = a.alu + ".plus.tmp"; minus_bed = a.alu + ".minus.tmp"
    with open(a.alu) as fh, open(plus_bed,"w") as pp, open(minus_bed,"w") as mm:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#","track","browser")): continue
            f = ln.rstrip("\n").split("\t")
            strand = f[5] if len(f) > 5 else "+"
            (pp if strand == "+" else mm).write("\t".join(f[:3])+"\n")

    # ref_cov and mismatch counts on the Alu SENSE strand
    ref_cov = {b:0 for b in "ACGT"}
    counts = {}  # 'X>Y' sense-strand directed mismatches

    def run(bed, minus):
        # samtools mpileup: -l regions, -Q baseq, -q mapq, -f ref. Output col5 = read bases.
        cmd = ["samtools","mpileup","-l",bed,"-f",a.fasta,
               "-Q",str(a.min_baseq),"-q",str(a.min_mapq),"-d","100000",a.bam]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        for line in p.stdout:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5: continue
            ref = c[2].upper()
            if ref not in "ACGT": continue
            bases = c[4]
            # count matches (.,) and each mismatch letter, ignoring indels/caps markers
            i=0; n=len(bases); acgt={"A":0,"C":0,"G":0,"T":0}; match=0
            while i < n:
                ch = bases[i]
                if ch in ".,": match+=1; i+=1
                elif ch=="^": i+=2  # skip mapping-qual char after ^
                elif ch=="$": i+=1
                elif ch in "+-":
                    j=i+1; num=""
                    while j<n and bases[j].isdigit(): num+=bases[j]; j+=1
                    i = j + (int(num) if num else 0)
                elif ch.upper() in "ACGT": acgt[ch.upper()]+=1; i+=1
                else: i+=1
            # genomic-forward ref base and observed counts -> flip to sense for '-' Alus
            comp={"A":"T","T":"A","G":"C","C":"G"}
            if not minus:
                sense_ref = ref
                obs = acgt
            else:
                sense_ref = comp[ref]
                obs = {comp[b]:acgt[b] for b in "ACGT"}
            # matches count as sense_ref
            ref_cov[sense_ref]+= match + obs[sense_ref]
            for b in "ACGT":
                if b != sense_ref and obs[b]>0:
                    counts[f"{sense_ref}>{b}"] = counts.get(f"{sense_ref}>{b}",0)+obs[b]
                    ref_cov[sense_ref]+=obs[b]
        p.wait()

    run(plus_bed, False)
    run(minus_bed, True)

    a_cov = ref_cov["A"]; ag = counts.get("A>G",0)
    aei = 100.0*ag/a_cov if a_cov else 0.0
    # noise floor: mean of A>C and A>T rates
    ac = counts.get("A>C",0); at = counts.get("A>T",0)
    noise = 100.0*(ac+at)/(2*a_cov) if a_cov else 0.0
    sn = aei/noise if noise>0 else float("inf")
    with open(a.out,"w") as o:
        o.write("sample\tAEI_percent\tAG_mismatches\tA_coverage\tsignal_to_noise\tnoise_floor_percent\n")
        o.write(f"{a.sample}\t{aei:.6f}\t{ag}\t{a_cov}\t{sn:.3f}\t{noise:.6f}\n")
    import os
    os.remove(plus_bed); os.remove(minus_bed)
    print(f"[AEI] {a.sample}: AEI={aei:.4f}% A>G={ag} A_cov={a_cov} S/N={sn:.2f}")

if __name__=="__main__":
    main()
