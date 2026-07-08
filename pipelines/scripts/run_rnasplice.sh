#!/usr/bin/env bash
# Run nf-core/rnasplice (rMATS, DEXSeq, edgeR, SUPPA2) on this Mac via conda profile.
# Responder-vs-nonresponder differential splicing on the melanoma ICB subset.
# Usage: run_rnasplice.sh <samplesheet.csv> <contrastsheet.csv> <outdir> [extra args]
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh"
# rnasplice 1.0.4's nextflow.config uses the legacy check_max()/max_memory pattern,
# rejected by Nextflow 26.x's v2 config parser — restore the v1 parser for this pipeline.
export NXF_SYNTAX_PARSER=v1
SS="${1:?samplesheet}"; CS="${2:?contrastsheet}"; OUT="${3:?outdir}"; shift 3 || true
REF="$REPO/reference/GRCh38"
WORK="$REPO/results/large/nf_work_rnasplice"
mkdir -p "$OUT" "$WORK"
ARGS=(
  -profile conda
  -c "$REPO/pipelines/conf/mac_arm64.config"
  -work-dir "$WORK"
  --input "$SS"
  --contrasts "$CS"
  --outdir "$OUT"
  --fasta "$REF/GRCh38.primary_assembly.genome.fa"
  --gtf "$REF/gencode.v46.primary_assembly.annotation.gtf"
  --aligner star
  --star_index "$REF/star_index"
  # rMATS is native arm64; DEXSeq/edgeR/SUPPA2 are R/python (arm64-fine)
  --rmats --dexseq_exon --edger_exon --suppa
  -resume
)
echo "nextflow run "$REPO/pipelines/nfcore/rnasplice-1.0.4/main.nf" ${ARGS[*]} $*"
nextflow run "$REPO/pipelines/nfcore/rnasplice-1.0.4/main.nf" "${ARGS[@]}" "$@"
