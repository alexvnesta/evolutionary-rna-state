#!/usr/bin/env python3
"""Build an nf-core/rnaseq samplesheet from a manifest (csv/json) chunk."""
import sys, csv, json
man, fq_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = json.load(open(man)) if man.endswith(".json") else list(csv.DictReader(open(man)))
with open(out, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample","fastq_1","fastq_2","strandedness"])
    for r in rows:
        run=r["run_accession"]; key=r["sample_title"].replace(" ","_")
        # strandedness = reverse (NOT auto). nf-core's `auto` triggers
        # FASTQ_SUBSAMPLE_FQ_SALMON, which needs fq=0.12.0 — a tool with NO
        # osx-arm64 conda build, so the env solve fails on this Mac. The pilot
        # (ERR2208952) empirically resolved to reverse: Salmon expected_format
        # ISR, RSeQC 92.3% antisense, strand-check status "pass". These are
        # TruSeq stranded libraries (Gide PRJEB23709 + Riaz PRJNA356761), all
        # reverse-stranded, so we set it explicitly and skip the arm64-broken
        # auto-detection subsample.
        w.writerow([key, f"{fq_dir}/{run}_1.fastq.gz", f"{fq_dir}/{run}_2.fastq.gz", "reverse"])
print(f"wrote {out} ({len(rows)} samples)")
