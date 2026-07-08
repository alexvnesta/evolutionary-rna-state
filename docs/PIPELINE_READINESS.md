# Pipeline Readiness — Differentiated Feature Modules

*What each differentiated/antigen module expects, so the moment a pipeline matrix
lands it runs without adaptation. All 6 entrypoints audited present + importable.*

`pilot_ingest.py` reads `results/features/` and runs every module whose input is
present, degrading gracefully (absent -> NA, never imputed). Drop these files into
`results/features/` and re-run the orchestrator.

## Expression-derived (computable NOW from gene TPM — already wired + validated)

| Module | Entrypoint | Input file | Status |
|---|---|---|---|
| GEP baseline | `pilot_gep.to_symbol_gene_matrix` + `gep_scores.score_all` | `quant_gene_tpm.parquet` | RUNS (n=52) |
| Regulator activity | `regulator_activity.build_regulator_activity` | `quant_gene_tpm.parquet` | RUNS (n=52) |

TPM layout accepted: genes×samples with a `gene_name` symbol column + Ensembl index
and run-accession sample columns (the 3-cohort pilot layout), OR samples×ENSG-columns.

## De-novo antigen (pending pipeline matrices — wired, awaiting input)

| Module | Entrypoint | Input file | Input schema |
|---|---|---|---|
| TE antigen | `differentiated.te_antigen.build_te_antigen_table(rows)` | `te_locus.parquet` | per-locus expressed-TE rows (family, locus, TPM, sequence/coords) + genome FASTA |
| Splicing neoantigen | `differentiated.splicing_neoantigen.build_feature_table(sample_junctions, sample_hla, fasta, exon_index)` | `splicing_junctions.parquet` | STAR SJ.out.tab junctions per sample + genome FASTA + exon index |
| Fusion antigen | `differentiated.fusion_antigen.build_fusion_feature_table(rows)` | `fusion_calls.parquet` | Arriba/STAR-Fusion calls (gene pair, breakpoints, frame) |
| Intron retention | `differentiated.intron_retention.compute_retained_intron_load(ir_data)` | `intron_retention.parquet` | per-intron retention ratios (IRFinder-S) |
| RNA editing (AEI) | `differentiated.rna_editing.compute_alu_editing_index(aei_data)` | `rna_editing_aei.parquet` | per-sample Alu editing index (REDItools) |
| SNV/indel neoantigen | `baseline.snv_indel_neoantigen.variants_from_maf(maf, proteome)` + `snv_indel_neoantigen_burden` | `variants.maf` | somatic MAF + proteome FASTA |

All antigen burdens feed the ONE shared MHC engine (`antigen_core.mhc_binding.count_binders`,
MHCflurry 2.2.0) fed by arcasHLA HLA typing (`hla_typing.parquet`). HLA table required
for antigen-burden modules; absent -> those features NA.

## Crosswalk (resolved)

Run-accession -> study sampleId -> response labels via `pilot_crosswalk`:
- Preferred: pipeline `run_catalog.csv` (all 3 cohorts incl. hugo) — 52/52 mapped.
- Fallback: ENA metadata (gide + riaz). Independently validated 40/40 agreement
  with the pipeline session's `combined_id_crosswalk.csv`.

## Test status
- `analysis/test_pilot_crosswalk.py` — 6/6 (crosswalk rules, catalog path, regulator scoring)
- `analysis/integration_test.py` — 11/11 (all module chains compose, cold-start)
