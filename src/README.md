# src — modeling core

The leakage-guarded modeling backbone for the evolutionary RNA-state project.
Consumes the harmonized analysis frame and de-novo phenotype matrices, constructs
the latent evolutionary RNA-state *S*, and organizes response prediction.

## Modules
- `data.py` — build the harmonized per-sample analysis frame across the five ICB cohorts.
- `features.py` — consume the de-novo multi-phenotype matrices (intron retention, RNA editing, TE/ERV, splicing, fusion) produced by the pipeline arm.
- `model.py` — latent evolutionary RNA-state *S* construction and the fold-contained response-prediction backbone (`latent_state`, `response_organization`).
