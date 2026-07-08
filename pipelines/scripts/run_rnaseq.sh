#!/usr/bin/env bash
# Run nf-core/rnaseq (HISAT2 genome BAMs + Salmon quant) on this Mac via conda.
# Usage: run_rnaseq.sh <samplesheet.csv> <outdir> [extra nextflow args...]
#
# ALIGNER = HISAT2, not STAR. The only osx-arm64 conda build of STAR 2.7.11b
# (the latest STAR release; builds haf7d672_6/_7/_8) is broken: it opens FASTQ
# input but ingests ZERO reads and writes empty BAMs. Reproduced exhaustively in
# this build environment — every arm64 build, the osx-64 build under Rosetta 2,
# stdin input, prebuilt + freshly-built + tiny self-built indexes, compressed and
# plain, single- and paired-end, all filters disabled — STAR always reads 0
# reads. HISAT2 2.2.1 (arm64-native conda) ingests the same reads correctly
# (verified: 100k real pilot reads read 100%). HISAT2 is splice-aware and emits
# the same coordinate-sorted genome BAMs the TE/ERV, intron-retention and
# RNA-editing subworkflows consume. Expression quantification still comes from
# Salmon pseudo-alignment (--pseudo_aligner salmon), which works on arm64
# (pilot: 42.2M reads, 85.1% mapped). Re-enable STAR on amd64/Docker.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh"
SS="${1:?samplesheet required}"; OUT="${2:?outdir required}"; shift 2 || true
REF="$REPO/reference/GRCh38"
FASTA="$REF/GRCh38.primary_assembly.genome.fa"
GTF="$REF/gencode.v46.primary_assembly.annotation.gtf"
WORK="$REPO/results/large/nf_work_rnaseq"
mkdir -p "$OUT" "$WORK"

ARGS=(
  -profile conda
  -c "$REPO/pipelines/conf/mac_arm64.config"
  -c "$REPO/pipelines/conf/arm64_module_overrides.config"
  -work-dir "$WORK"
  --input "$SS"
  --outdir "$OUT"
  --fasta "$FASTA"
  --gtf "$GTF"
  --aligner hisat2             # arm64-native; produces genome BAMs for subworkflows
  --pseudo_aligner salmon      # expression quant (STAR-free, works on arm64)
  # save_reference + skip_bbsplit set (typed) in mac_arm64.config, not here
  -resume
)
# HISAT2 index: nf-core builds it. The module's HISAT2_BUILD auto-degrades to a
# PLAIN index (no splice sites / exons) when available RAM < params.hisat2_build_memory
# (default 200 GB). On this 64 GB Mac it therefore builds the plain graph index
# (~5-6 GB RAM), which still yields valid coordinate-sorted genome BAMs. Persisted
# via --save_reference (set in mac_arm64.config) and reused on -resume.
echo "nextflow run "$REPO/pipelines/nfcore/rnaseq-3.26.0/main.nf" ${ARGS[*]} $*"
nextflow run "$REPO/pipelines/nfcore/rnaseq-3.26.0/main.nf" "${ARGS[@]}" "$@"
