#!/usr/bin/env bash
# Standalone, mount-resilient Alu Editing Index (AEI) runner for the cohort.
#
# WHY THIS EXISTS: the Nextflow rna_editing run repeatedly had its head process
# killed by the intermittent host-mount drop during the AEI phase (heavy BAM I/O
# correlates with the drop). AEI is per-sample independent, so a plain parallel
# bash loop with per-sample retry is far more robust: a mount blip only costs one
# sample a retry instead of killing the whole orchestration. Each compute_aei.py
# call is the exact command the ALU_EDITING_INDEX Nextflow process would run.
#
# Usage: run_aei_standalone.sh <samplesheet.csv sample,bam,bai> <outdir> [NPAR]
set -uo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SS="${1:?samplesheet (sample,bam,bai)}"
OUT="${2:-$REPO/results/rna_editing_cohort/rna_editing/aei}"
NPAR="${3:-6}"                      # concurrent samples (AEI is single-threaded each)
EENV="$REPO/.nextflow_home/conda_cache/env-dba0831ccec915cf4db49bd8d700adf7"
PY="$EENV/bin/python"
BIN="$REPO/pipelines/rna_editing/bin/compute_aei.py"
FASTA="$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa"
# AEI is a POOLED genome-wide ratio (sum A>G / sum A-coverage over Alus). Full 1.18M-Alu
# pileup is ~5.4 h/sample (27 h for 40 even 8-way) because pysam does a per-interval
# pileup (~16 ms each). A FIXED representative 200k-Alu subset (seed=42, all 24 chroms,
# ~millions of A positions/sample) applied identically to all samples gives fully
# comparable cross-sample AEI at ~55 min/sample (~3 h for 40 at 12-way). Set
# AEI_FULL=1 to use the complete Alu set instead.
if [ "${AEI_FULL:-0}" = "1" ]; then
  ALU="$REPO/results/large/alu.hg38.bed6"
else
  ALU="$REPO/results/large/alu.hg38.sub200k.bed6"
fi
MINQ=60; MINBQ=25                   # HISAT2 unique-mapper MAPQ=60; baseq 25 (matches editing.config)
mkdir -p "$OUT"

[ -x "$PY" ]      || { echo "[error] editing env python missing: $PY"; exit 1; }
[ -f "$BIN" ]     || { echo "[error] compute_aei.py missing: $BIN"; exit 1; }
[ -f "$ALU" ]     || { echo "[error] Alu BED missing: $ALU"; exit 1; }
[ -f "$FASTA.fai" ] || "$EENV/bin/samtools" faidx "$FASTA"

# one sample, with retry (mount drops mid-read -> retry the sample)
run_one() {
  local sample="$1" bam="$2"
  local out="$OUT/${sample}.aei.tsv"
  # skip if already done and non-empty (idempotent resume)
  if [ -s "$out" ] && [ "$(wc -l < "$out")" -ge 2 ]; then
    echo "[skip] $sample (already done)"; return 0
  fi
  local tries=0 max=4
  while [ $tries -lt $max ]; do
    tries=$((tries+1))
    if "$PY" "$BIN" --bam "$bam" --fasta "$FASTA" --alu "$ALU" \
         --sample "$sample" --min-baseq "$MINBQ" --min-mapq "$MINQ" \
         --out "$out.tmp" 2> "$OUT/${sample}.log"; then
      mv "$out.tmp" "$out"
      echo "[ok  ] $sample ($(tail -1 "$out" | cut -f2) % AEI)"
      return 0
    fi
    echo "[retry] $sample attempt $tries/$max failed; backoff"
    rm -f "$out.tmp"; sleep $((tries*15))
  done
  echo "[FAIL] $sample after $max attempts"; return 1
}
# drive N samples in parallel — read the whole samplesheet into arrays first
# (avoid piping into the loop: a pipe subshell + background jobs is fragile)
echo "[aei] $(date '+%H:%M:%S') start, NPAR=$NPAR, out=$OUT"
SAMPLES=(); BAMS=()
while IFS=, read -r sample bam bai _; do
  sample="${sample%$'\r'}"; bam="${bam%$'\r'}"
  [ "$sample" = "sample" ] && continue    # header
  [ -z "$sample" ] && continue
  SAMPLES+=("$sample"); BAMS+=("$bam")
done < "$SS"
echo "[aei] ${#SAMPLES[@]} samples to process"

for i in "${!SAMPLES[@]}"; do
  # throttle to NPAR concurrent jobs
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 5; done
  run_one "${SAMPLES[$i]}" "${BAMS[$i]}" &
done
wait
echo "[aei] $(date '+%H:%M:%S') all samples done"

# merge (only if at least one per-sample TSV exists)
shopt -s nullglob
TSVS=("$OUT"/*.aei.tsv)
shopt -u nullglob
NDONE=${#TSVS[@]}
echo "[aei] $NDONE per-sample AEI produced; merging"
if [ "$NDONE" -gt 0 ]; then
  "$PY" "$REPO/pipelines/rna_editing/bin/merge_aei.py" \
    --out "$OUT/cohort_aei.tsv" "${TSVS[@]}" && \
    echo "[aei] cohort_aei.tsv written: $(wc -l < "$OUT/cohort_aei.tsv") lines"
else
  echo "[aei] no per-sample TSVs — nothing to merge"
fi
