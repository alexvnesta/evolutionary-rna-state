# Data inventory

All datasets below are **public**. This repo commits metadata only (manifests,
receipts, run catalog, clinical labels); raw reads are downloaded to
git-ignored paths at analysis time.

## Cohorts

| Dataset | Cancer | N (RNA-seq) | Raw FASTQ | Response labels | Accession | Role |
|---------|--------|-------------|-----------|-----------------|-----------|------|
| **Gide 2019** | Melanoma | 91 (73 PRE) | ENA (open) | RECIST + arm/timepoint | `PRJEB23709` | Development |
| **Riaz 2017** | Melanoma | 109 (Pre+On) | SRA/ENA (open) | RECIST, GSM-mapped | `PRJNA356761` / `GSE91061` | External validation |
| Zhao/Cloughesy 2019 | Glioblastoma | 34 (+57 WXS) | ENA (open) | anti-PD-1 response | `PRJNA482620` | Extension |
| Kim 2018 | Gastric | 78 (+110 WXS) | ENA (open) | pembrolizumab response | `PRJEB25780` | Extension |
| IMvigor210 | Urothelial | ~348 (matrix) | **EGA controlled** (raw) | atezo RECIST/OS | processed R pkg | Extension (processed only) |
| TISCH2 | Pan-cancer | millions of cells | N/A (processed scRNA) | many IO cohorts | tisch.comp-genomics.org | scRNA validation |

Sizes and preprocessing notes are in `data/dataset_summary.csv`.

## Primary axis: melanoma ICB (raw reads available)

Gide 2019 and Riaz 2017 are the two melanoma cohorts with **open raw reads** and
RECIST response, making them the development/validation pair for the raw-read
components of the thesis. Pretreatment (PRE) samples are the primary analysis
set; on-treatment samples are inventoried for longitudinal extension.

## Committed metadata files

- `data/dataset_summary.csv` — the table above, machine-readable.
- `data/catalog/run_catalog.csv` — 200 sequencing runs with accessions,
  patient/arm/timepoint, library metadata, FASTQ URLs + MD5s, and joined
  clinical fields (RECIST, responder/non-responder, therapy, OS, vital).
- `data/manifests/selection_manifest.csv` — runs selected for analysis, tagged
  by cohort role.
- `data/manifests/pilot_manifest.json` — pilot subset for pipeline bring-up.
- `data/manifests/download_receipts.json` — per-file download receipts with
  observed vs expected MD5 and byte counts (integrity provenance).
- `data/response/Melanoma-GSE91061.clinical.tsv` — Riaz 2017 clinical/response.
- `data/response/Melanoma-PRJEB23709.clinical.tsv` — Gide 2019 clinical/response.

## Access & ethics

Open datasets only in this repo. Controlled-access raw data (e.g. IMvigor210 via
EGA) is never committed; only its publicly distributed processed matrix is used.
Raw reads and derived large files are git-ignored.
