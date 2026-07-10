#!/bin/bash
# arcasHLA Class-I typing over accession-named cohort BAMs (presentation-layer step 1).
# Reproduces the known-good ERR2208909 invocation. No GPU. Class I only (A,B,C) for NetMHCpan.
set -uo pipefail
ENV=/Users/alex/.claude-science/conda/envs/antigen-hla-test/bin
K046=/Users/alex/.claude-science/conda/envs/kallisto046/bin   # arcasHLA 0.5.0 needs kallisto 0.46.x (pseudoalignments.tsv); 0.52 writes bus format instead
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
SCRIPTS="$REPO/tools/arcasHLA/scripts"   # vendored scripts (conda pkg has broken import; wrapper uses blocked realpath)
export PATH="$K046:$ENV:$PATH"           # kallisto046 FIRST so `kallisto` resolves to 0.46.1
export PYTHONPATH="$SCRIPTS:${PYTHONPATH:-}"
PY="$ENV/python3"                        # antigen-hla-test python (has numpy/scipy); kallisto046 python does not
BAMDIR="$REPO/results/editing_bams"
OUT="$REPO/results/hla"
mkdir -p "$OUT"
cd "$REPO"

# singleton lock — refuse to run if another copy is active (mkdir is atomic)
LOCK="$OUT/.hla_typing.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "ANOTHER HLA BATCH IS RUNNING (lock $LOCK exists) — exiting"; exit 3
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

BAMS=$(ls "$BAMDIR"/*.hisat2.sorted.bam)
N=$(echo "$BAMS" | wc -l | tr -d ' ')
i=0
for bam in $BAMS; do
  i=$((i+1))
  acc=$(basename "$bam" | sed 's/\.hisat2.*//')
  od="$OUT/$acc"
  if [ -f "$od/${acc}.genotype.json" ]; then echo "[$i/$N] $acc SKIP (done)"; continue; fi
  mkdir -p "$od"
  echo "[$i/$N] $acc START $(date +%H:%M:%S)"
  [ -f "${bam}.bai" ] || samtools index -@ 8 "$bam"
  # 1) extract chr6 HLA reads -> paired fastqs (direct script call); reuse if already extracted
  fq1=$(ls "$od"/*.extracted.1.fq.gz 2>/dev/null | head -1)
  fq2=$(ls "$od"/*.extracted.2.fq.gz 2>/dev/null | head -1)
  if [ -z "$fq1" ] || [ -z "$fq2" ]; then
    "$PY" "$SCRIPTS/extract.py" "$bam" -o "$od" -t 8 -v > "$od/${acc}.extract.log" 2>&1
    fq1=$(ls "$od"/*.extracted.1.fq.gz 2>/dev/null | head -1)
    fq2=$(ls "$od"/*.extracted.2.fq.gz 2>/dev/null | head -1)
  fi
  if [ -z "$fq1" ] || [ -z "$fq2" ]; then echo "[$i/$N] $acc FAIL extract"; continue; fi
  # 2) genotype Class I (A,B,C) — direct script call
  "$PY" "$SCRIPTS/genotype.py" "$fq1" "$fq2" -g A,B,C -o "$od" -t 8 -v > "$od/${acc}.genotype.log" 2>&1
  # arcasHLA names genotype json off the fq1 stem; normalise to <acc>.genotype.json
  gj=$(ls "$od"/*.genotype.json 2>/dev/null | head -1)
  if [ -n "$gj" ]; then
    [ "$gj" = "$od/${acc}.genotype.json" ] || cp "$gj" "$od/${acc}.genotype.json"
    echo "[$i/$N] $acc DONE $(cat "$od/${acc}.genotype.json")"
  else
    echo "[$i/$N] $acc FAIL genotype"
  fi
done
echo "ALL DONE $(date +%H:%M:%S)"
