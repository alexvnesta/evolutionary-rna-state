# analysis/differentiated/bin — per-phenotype CLIs

Command-line drivers called by the differentiated (per-RNA-phenotype) Nextflow modules.

## Modules
- `ir_antigen_features.py` — turn the pipeline's intron-retention output into antigen features.
- `build_te_antigen_table.py` — aggregate per-sample TE/ERV antigen burden into a cohort table.
- `fusion_burden_cli.py` — per-sample fusion-neoantigen-burden driver.
