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

## Derived data joins (for the analysis sessions)

Two joins that unblock analysis directions named in `HANDOFF_rna_state_next.md`.
Provenance and method detail: `_joins_provenance.json`.

- **`riaz_clonality.csv`** — per-patient clonal architecture reconstructed from
  the Riaz Table S3 somatic MAF (tumor VAF `Taf`). For each of 68 patients:
  clonal VAF peak (2-component GMM, or 90th-percentile fallback at <20
  mutations), `n_clonal`/`n_subclonal` (CCF≥0.8 threshold), `subclonal_fraction`,
  continuous `mean_ccf`, and a `heterogeneity_index` (CCF spread, an
  intratumoral-heterogeneity proxy). Joins to the Riaz RNA cohort on the shared
  `Pt##` patient namespace: **51 patients carry both RNA-seq and a labeled
  response (11 R / 40 N)** — the test set for the RNA-by-clonality interaction.
  Caveat: no published FACETS purity in the supplements, so CCF is
  purity-uncorrected; the continuous metrics are primary and the clonal call is
  coarse. Recompute with PyClone-VI + purity for production.
- **`gide_id_crosswalk.csv`** — run-level map from ENA PRJEB23709 identifiers
  (`run_accession`, `sample_title` = `ipiPD1_N`/`PD1_N`) to iAtlas
  `iatlas_sampleId` (`PD0N_Pre`, `iPiN_On`). Lets expression keyed on
  SRA/trial IDs join to the cBioPortal burden/neoantigen features. Built as a
  deterministic (arm, patient-number) map and **validated by clinical
  fingerprint: 75/75 patients concordant on treatment arm, timepoint structure,
  and response; all 91 runs round-trip to valid iAtlas sample IDs.**
- **`hugo_clonality.csv`** — per-patient clonal architecture for Hugo 2016 (38
  patients, 21 R / 17 N), computed as a **purity-corrected CCF**:
  `CCF = clip(VAF·2/purity, 0, 1)` with VAF parsed from the S1D MAF caller
  fields (validated r=0.996 against VarScan2's stated %) and purity from the
  published S1B table. Same schema as `riaz_clonality.csv`. Mutation counts
  track the published `TotalNonSyn` (r=0.947). Caveats: 10 patients have
  purity<0.30 (`low_purity_flag`), where the correction saturates and over-calls
  clonal; multiplicity is fixed at 1 with no allele-specific CNA correction.
  **Not cross-cohort comparable to `riaz_clonality.csv` on absolute clonality** —
  Riaz uses a purity-free VAF/clonal-peak proxy (median subclonal fraction 0.70),
  Hugo uses purity-corrected CCF (0.25); the difference is method, not biology.
  Use each within-cohort only. Joins to the Hugo RNA cohort
  (GSE78220/PRJNA356839) on the `Pt##` namespace once that cohort is catalogued.
