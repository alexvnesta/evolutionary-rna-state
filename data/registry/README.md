# Dataset registry

Organized record of the datasets, choices, and supplemental data behind this
project — and the source for the manuscript's **"Data — what we used and why"**
section. All numbers here are verified against source (ENA run catalog +
cBioPortal live pull), not estimated.

## Files

| File | Contents |
|------|----------|
| `DATASET_REGISTRY.xlsx` | Formatted workbook, 5 sheets + README (open this for the human-readable master). |
| `cohorts.csv` | Master table — 5 melanoma ICB cohorts, roles, RNA/WES availability, rationale. |
| `rna_phenotypes.csv` | The 6 thesis RNA phenotypes → data source per phenotype. |
| `covariates.csv` | Confounders-to-beat (TMB, purity, ploidy, clonality, composition) + benchmark signatures. |
| `supplements.csv` | Per-paper supplemental-data inventory + access status. |
| `coverage_verified.csv` | Per-cohort field coverage, pulled live from cBioPortal. |

Per-sample clinical/genomic values pulled from cBioPortal live are in
`../cbioportal/*.clinical.tsv` (one file per cohort).

## The cohort set at a glance

**RNA raw reads (FASTQ, from ENA)** — the thesis core, because annotation-based
pipelines discard the non-reference signal we are after:

- **Gide 2019** (`PRJEB23709`) — RNA development lead. Deepest/longest reads
  (median 48.5M pairs / 9.9 Gbp, HiSeq 2500) and balanced response
  (CR17/PR32/SD13/PD29). No WES.
- **Riaz 2017** (`PRJNA356761`) — clonal-evolution + multimodal anchor. The only
  cohort with longitudinal WES + RNA + TCR and published PyClone CCF / clonality.
- **Hugo 2016** (`PRJNA356839` / `GSE78220`) — third RNA cohort. Both modalities;
  source of the IPRES resistance signature we benchmark against.

**Variant calls only (cBioPortal iAtlas-harmonized MAF; no raw exomes needed):**

- **Riaz + Hugo** — matched to the RNA cohorts (true TMB, driver mutations).
- **Liu 2019** (`mel_iatlas_liu_2019`, 122 samples) — largest open melanoma ICB
  variant set, for driver/TMB statistical power.
- **DFCI 2019** (`mel_dfci_2019`, 144 samples) — the only cohort carrying
  clonal/subclonal mutation counts + purity + ploidy + GISTIC CNA.

## Why this set (short version)

1. **RNA raw reads + same tissue.** Only melanoma provides ≥3 open-raw-read ICB
   cohorts, so external validation is not confounded by tissue-of-origin.
2. **Split leads by strength.** Gide has the best RNA substrate for de-novo
   phenotype quantification; Riaz has the deep multimodal/clonal data; Hugo
   shores up the (otherwise responder-thin) validation side.
3. **Variant side without raw exomes.** cBioPortal's iAtlas harmonization gives
   one-pipeline somatic calls + TMB across cohorts, so the driver/clonal arm is
   backed by published calls with no dbGaP/EGA access wall and no cross-cohort
   caller batch effect.
4. **A bonus alignment.** iAtlas already quantifies neoantigens by mechanism —
   `SPLICE_`, `ERV_`, `FUSION_` — i.e. three of our six RNA phenotypes,
   per-sample, at ~99–100% coverage for Gide/Hugo/Liu (68% for Riaz). This
   gives an independent, orthogonal test of the **co-variation claim** before
   our own raw-read pipeline is stood up.

## Verified coverage highlights (live cBioPortal pull)

- Splice/ERV/Fusion neoantigen categories: **Gide 99%, Hugo 100%, Liu 100%,
  Riaz 68%** (DFCI: not provided — carries clonal/CNA instead).
- True TMB (`TMB_NONSYNONYMOUS`): Riaz/Hugo/Liu/DFCI. Purity+ploidy+clonal
  counts: DFCI only.
- Response (RECIST): all five (Gide/Hugo/Liu 100%, Riaz 92%).

See `coverage_verified.csv` for the full matrix.
