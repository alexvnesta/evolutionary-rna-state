# HLA-I typing — presentation-layer step 1 (Phase 1)

**Date:** 2026-07-10 | **Status:** COMPLETE (16/16 accession-named cohort BAMs) | **Compute:** local, no GPU.

## What ran
arcasHLA 0.5.0 `extract` → `genotype -g A,B,C` on the 16 accession-named BAMs in `results/editing_bams/`
(10 Gide / 2 Hugo / 4 Riaz — the exact set the non-reference matrix covers). Output per sample in
`results/hla/<acc>/<acc>.genotype.json`; aggregated to `results/predictor/hla_alleles.parquet` (+ .csv).

## Environment gotchas (reproducibility — these cost real debugging time)
1. **arcasHLA entry point:** the conda `arcasHLA` (env `antigen-hla-test`) has a broken `arcas_utilities`
   import; the vendored wrapper `tools/arcasHLA/arcasHLA` calls `realpath` (sandbox-blocked). FIX: call the
   vendored scripts directly — `python3 tools/arcasHLA/scripts/{extract,genotype}.py` with
   `PYTHONPATH=tools/arcasHLA/scripts`.
2. **kallisto version:** arcasHLA 0.5.0 needs `kallisto pseudo` to write `pseudoalignments.tsv`. kallisto
   0.52 (default on PATH) DEPRECATED `pseudo` and writes BUS format instead → genotype crashes with
   `FileNotFoundError: pseudoalignments.tsv`. FIX: put env `kallisto046` (kallisto 0.46.1) FIRST on PATH.
3. **python for the scripts:** kallisto046's python lacks numpy; use `antigen-hla-test/bin/python3` to run the
   scripts while `kallisto` resolves to 0.46.1. Runner `analysis/run_hla_typing.sh` encodes all three fixes
   and a singleton lock.

## Result summary
- 16/16 typed. Mean HLA-I heterozygosity fraction (distinct 2-field alleles / 6) = **0.885**; 9/16 fully
  heterozygous at all three loci. Heterozygosity is stored as the Chowell-2018 covariate `hla_het_fraction`.
- 59 distinct 2-field alleles across the cohort.

## KNOWN LIMITATION for the binding step (must handle before NetMHCpan/MHCflurry)
**9/16 samples carry >=1 allele NOT directly in mhcflurry's supported list** (`unsupported` column). Two
causes: (a) high-resolution / rare 2-field codes (e.g. A*66:57, B*40:538) not in the mass-spec training set;
(b) **null alleles** (suffix N/Q, e.g. B*44:345N, A*02:945N, B*13:07N) — these are NOT EXPRESSED proteins and
must be DROPPED from presentation (a null allele presents nothing), not mapped to a nearest neighbour. The
binding step must: drop null/low-expression (N/Q/L) alleles, and for genuinely-unsupported expressed alleles
either use NetMHCpan (broader allele coverage) or a documented nearest-supported fallback — logged per sample.
This is a real design item, not a typing failure; the calls themselves are valid arcasHLA output.
