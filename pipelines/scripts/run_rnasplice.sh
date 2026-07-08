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

# SOURCE selects the input mode:
#   genome_bam (DEFAULT on this arm64 Mac) — feed pre-aligned HISAT2 spine BAMs,
#     bypassing the broken arm64 STAR entirely. Enables rMATS + DEXSeq + edgeR.
#     SUPPA needs salmon_results and is NOT available from genome_bam, so it is
#     only added when SOURCE=fastq.
#   fastq — the upstream default (STAR-align then split); use on amd64/Docker.
SOURCE="${SOURCE:-genome_bam}"

ARGS=(
  -profile conda
  -c "$REPO/pipelines/conf/mac_arm64.config"
  -c "$REPO/pipelines/conf/rnasplice_arm64_overrides.config"
  -work-dir "$WORK"
  --input "$SS"
  --contrasts "$CS"
  --outdir "$OUT"
  --source "$SOURCE"
  --fasta "$REF/GRCh38.primary_assembly.genome.fa"
  --gtf "$REF/gencode.v46.primary_assembly.annotation.gtf"
)
if [ "$SOURCE" = "fastq" ]; then
  # full path: STAR align (amd64/Docker only — arm64 STAR ingests 0 reads) + all tools
  ARGS+=( --aligner star --star_index "$REF/star_index"
          --rmats --dexseq_exon --edger_exon --suppa )
else
  # genome_bam path: no alignment; splice tools that consume genome BAM directly.
  # DEXSeq + edgeR are R-based and run cleanly on arm64 conda.
  # SUPPA needs salmon_results (fastq source) — not available here.
  # MISO sashimi plots are disabled: misopy=0.5.4 is python-2.7-only with no
  # osx-arm64 build, and the plots are optional visualisation, not results.
  # rMATS is OPT-IN (RMATS=1): its env (rmats=4.3.0 + r-pairadise) has an
  # unsatisfiable R-toolchain conflict on osx-arm64 (r-base 4.2 vs 4.5, gsl,
  # libgfortran). Off by default; enable only on amd64/Docker. See RUNBOOK §5.
  # rnasplice defaults rmats/suppa/dexseq/edger all to true in nextflow.config,
  # so unwanted tools must be EXPLICITLY disabled with `<tool> false`.
  ARGS+=( --dexseq_exon --edger_exon --suppa false --sashimi_plot false )
  if [ "${RMATS:-0}" = "1" ]; then ARGS+=( --rmats ); else ARGS+=( --rmats false ); fi
fi
ARGS+=( -resume )
echo "nextflow run "$REPO/pipelines/nfcore/rnasplice-1.0.4/main.nf" ${ARGS[*]} $*"
nextflow run "$REPO/pipelines/nfcore/rnasplice-1.0.4/main.nf" "${ARGS[@]}" "$@"
