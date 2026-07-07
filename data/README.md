# Data pointers

**No raw data lives here.** This directory holds *manifests* describing where data
lives and how to obtain it.

For each dataset add a file `data/<dataset_id>.manifest` (or a subdirectory with a
`README.md`) recording:

- Source / repository (e.g. GEO, SRA, ENA, dbGaP, TCGA)
- Accession(s)
- Access level (open vs controlled) and any DUA/approval required
- Assay(s): RNA-seq, WES/WGS, etc.
- Relevant clinical/outcome fields (e.g. ICB response, survival)
- Checksums where practical

Controlled-access data must never be committed, even in part.
