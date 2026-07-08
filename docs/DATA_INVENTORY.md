# Data inventory

All datasets below are **public**. This repo commits metadata only (manifests,
receipts, run catalog, clinical labels, harmonized cBioPortal pulls); raw reads
are downloaded to git-ignored paths at analysis time.

**The authoritative, machine-readable registry is in `data/registry/`**
(`DATASET_REGISTRY.xlsx` + per-sheet CSVs). This page is the prose summary.

## Finalized cohort set (melanoma ICB)

### RNA raw reads (FASTQ, from ENA) — the thesis core

Annotation-based pipelines discard non-reference signal (intron retention, TE
transcripts, editing, novel junctions), which is exactly where the evolutionary
RNA-state fingerprint lives — so the core requires raw reads. Melanoma is the
only tumor type with ≥3 open-raw-read ICB cohorts, keeping external validation
free of tissue-of-origin confounding.

| Cohort | RNA runs | Depth (median) | Instrument | WES | cBioPortal | Role |
|--------|----------|----------------|------------|-----|------------|------|
| **Gide 2019** `PRJEB23709` | 91 (73 PRE) | 48.5M / 9.9 Gbp | HiSeq 2500 | none | `mel_iatlas_gide_2019` | RNA development lead |
| **Riaz 2017** `PRJNA356761` | 109 (51 PRE) | 52.3M / 5.3 Gbp | Genome Analyzer | 150× + PyClone CCF | `mel_iatlas_riaz_nivolumab_2017` | Clonal + multimodal anchor |
| **Hugo 2016** `PRJNA356839`/`GSE78220` | 24 (PRE) | (GEO) | HiSeq | 140× | `mel_iatlas_hugo_ucla_2016` | Third RNA cohort (validation) |

Split-lead rationale: Gide has the deepest/longest reads and balanced response
(best substrate for de-novo phenotype quantification and feature development);
Riaz has the only longitudinal WES+RNA+TCR data and published clonal/CCF
analysis (clonal-evolution arm + true-TMB confounder); Hugo adds a second
validation cohort and the IPRES resistance signature.

### Variant calls only (cBioPortal iAtlas MAF — no raw exomes)

Per the scope decision, the driver/clonal arm uses **published somatic calls**,
not self-called raw exomes. cBioPortal's iAtlas harmonization applies one
pipeline across cohorts (no cross-cohort caller batch effect) and avoids the
dbGaP/EGA access wall on raw WXS.

| Cohort | Samples | Carries | cBioPortal |
|--------|---------|---------|------------|
| Riaz 2017 | 73 (MAF) | TMB, drivers, neoantigen categories | `mel_iatlas_riaz_nivolumab_2017` |
| Hugo 2016 | 27 (MAF) | TMB, drivers, neoantigen categories | `mel_iatlas_hugo_ucla_2016` |
| **Liu 2019** | 122 (MAF) | TMB, drivers (largest open set) | `mel_iatlas_liu_2019` |
| **DFCI 2019** | 144 (MAF+CNA) | clonal/subclonal counts, purity, ploidy, GISTIC | `mel_dfci_2019` |

### Bonus: mechanistic neoantigen categories already harmonized

iAtlas quantifies neoantigens by mechanism — `SPLICE_NEOANTIGEN`,
`ERV_NEOANTIGEN`, `FUSION_NEOANTIGEN` — i.e. three of our six RNA phenotypes,
per-sample (Gide 99%, Hugo/Liu 100%, Riaz 68%). This gives an independent test
of the co-variation claim before the raw-read pipeline exists. Pulled to
`data/cbioportal/*.clinical.tsv`.

## Datasets considered and excluded

- **Van Allen 2015** — raw WXS access-controlled (dbGaP); usable only as
  variant calls, and the melanoma trio + Liu/DFCI already cover the variant side.
- **IMvigor210 (urothelial)** — raw = EGA controlled; processed counts only;
  different tissue.
- **Zhao/GBM, Kim/gastric** — open raw reads but different tissues; deferred to a
  cross-tissue extension (scope decision: melanoma depth first).
- **TISCH2** — processed scRNA; possible orthogonal validation, not core.

## Committed metadata files

- `data/dataset_summary.csv` — the table above, machine-readable.
- `data/catalog/run_catalog.csv` — 228 RNA-Seq runs (riaz2017 109, gide2019 91,
  hugo2016 28) with accessions, patient/arm/timepoint, library metadata, FASTQ
  URLs + MD5s, and joined clinical fields (RECIST, responder/non-responder,
  therapy, OS, vital). Hugo RNA cohort = SRP070710 / PRJNA312948 (GSE78220),
  cataloged from ENA with clinical fields from Hugo S1A.
- `data/manifests/selection_manifest.csv` — runs selected for analysis, tagged
  by cohort role.
- `data/manifests/pilot_manifest.json` — pilot subset for pipeline bring-up.
- `data/manifests/download_receipts.json` — per-file download receipts with
  observed vs expected MD5 and byte counts (integrity provenance).
- `data/response/Melanoma-GSE91061.clinical.tsv` — Riaz 2017 clinical/response.
- `data/response/Melanoma-PRJEB23709.clinical.tsv` — Gide 2019 clinical/response.
- `data/registry/` — **dataset registry** (`DATASET_REGISTRY.xlsx` + CSVs): the
  authoritative record of cohorts, phenotype→source map, covariates, supplement
  inventory, and verified coverage. Source for the manuscript Data section.
- `data/cbioportal/*.clinical.tsv` — per-sample harmonized clinical/genomic
  values pulled live from cBioPortal (TMB, neoantigen categories by mechanism,
  response, survival, purity/ploidy/clonality where available).

## Access & ethics

Open datasets only in this repo. Controlled-access raw data (e.g. IMvigor210 via
EGA) is never committed; only its publicly distributed processed matrix is used.
Raw reads and derived large files are git-ignored.
