#!/usr/bin/env python3
"""Build an nf-core/rnaseq samplesheet from a manifest (csv/json) chunk."""
import sys, csv, json
man, fq_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = json.load(open(man)) if man.endswith(".json") else list(csv.DictReader(open(man)))
with open(out, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample","fastq_1","fastq_2","strandedness"])
    for r in rows:
        run=r["run_accession"]; key=r["sample_title"].replace(" ","_")
        w.writerow([key, f"{fq_dir}/{run}_1.fastq.gz", f"{fq_dir}/{run}_2.fastq.gz", "auto"])
print(f"wrote {out} ({len(rows)} samples)")
