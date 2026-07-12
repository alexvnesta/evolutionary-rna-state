#!/bin/bash
set -eo pipefail
SHIM=/tmp/arcas_shim; mkdir -p $SHIM
ln -sf /Users/alex/.claude-science/conda/envs/kallisto046/bin/kallisto $SHIM/kallisto
ln -sf /Users/alex/.claude-science/conda/envs/antigen-hla-test/bin/python3 $SHIM/python3
ln -sf /Users/alex/.claude-science/conda/envs/nonref-callers/bin/samtools $SHIM/samtools
ln -sf /Users/alex/.claude-science/conda/envs/antigen-hla-test/bin/bedtools $SHIM/bedtools
ln -sf /Users/alex/.claude-science/conda/envs/antigen-hla-test/bin/pigz $SHIM/pigz
export PATH=$SHIM:$PATH
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
ARCAS=$REPO/tools/arcasHLA; export PYTHONPATH=$ARCAS/scripts
OUT=$REPO/results/hla_typing/gide_arcas; mkdir -p $OUT
cd $OUT
for TB in $REPO/results/rnaseq_cohort/hisat2/*PD1*.markdup.sorted.bam; do
  base=$(basename "$TB" .markdup.sorted.bam)
  if [ -f "${base}.markdup.sorted.extracted.genotype.json" ]; then echo "SKIP $base (done)"; continue; fi
  echo "=== $base $(date +%H:%M:%S) ==="
  T=/tmp/arcas_${base}; mkdir -p $T
  $ARCAS/arcasHLA extract "$TB" -o $OUT -t 6 --temp $T >/dev/null 2>&1 || { echo "EXTRACT FAIL $base"; continue; }
  $ARCAS/arcasHLA genotype ${base}.markdup.sorted.extracted.1.fq.gz ${base}.markdup.sorted.extracted.2.fq.gz \
     -g A,B,C -o $OUT -t 6 --temp $T >/dev/null 2>&1 || { echo "GENOTYPE FAIL $base"; continue; }
  rm -f ${base}.markdup.sorted.extracted.[12].fq.gz; rm -rf $T
  echo "DONE $base: $(cat ${base}.markdup.sorted.extracted.genotype.json 2>/dev/null)"
done
echo "=== ALL DONE $(date +%H:%M:%S) ==="
ls *.genotype.json | wc -l
