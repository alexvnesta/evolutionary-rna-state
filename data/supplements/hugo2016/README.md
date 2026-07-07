# Hugo et al. 2016 — supplemental data

**Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in
Metastatic Melanoma.** Hugo W, et al. *Cell* 165(1):35–44, 2016.
PMID 26997480. RNA-seq: GEO GSE78220. Anti-PD-1 (pembrolizumab/nivolumab).

Source files (publisher supplemental): `mmc1.xls` (Table S1 workbook),
`mmc2.xlsx` (Table S2 workbook), `mmc3.pdf` (supplemental figures).
MD5s in `raw/checksums.json`.

## Cleaned tables (`clean/`)

| CSV | Source | Rows×Cols | Contents |
|-----|--------|-----------|----------|
| `hugo2016_S1A.csv` | mmc1 S1A | 42×22 | Patient/sample characteristics: irRECIST, site, gender, age, OS, vital, prior MAPKi, anatomy, BRAF/NRAS status. |
| `hugo2016_S1B.csv` | mmc1 S1B | 39×45 | **WES + neoepitope calls per sample**: BRCA2-mutant flag, IPRES signature, Response, Purity/Ploidy (Sequenza), AvgCov, TotalNonSyn, indel load, MHC binder counts, mutation spectrum (A>G, C>T, …). |
| `hugo2016_S1C.csv` | mmc1 S1C | 37×6 | Neoepitope peptides: sequence, HLA type, gene, AA mutation, position. |
| `hugo2016_S1D.csv` | mmc1 S1D | **25393×17** | **Full somatic mutation table (MAF)**: Chr, Pos, NucMut, Sample, Gene, Transcript, MutType, AAmut, codon change. 38 samples, 9,877 genes. |
| `hugo2016_S1E.csv` | mmc1 S1E | 36×11 | Genes differentially mutated R vs NR: ratio, hit counts, p/adj, log-odds. |
| `hugo2016_S2A.csv` | mmc2 S2A | 693×7 | Differentially expressed genes (R vs NR): Mann-Whitney p, FDR, U-stat, diffAvg, avg.R, avg.NR. |
| `hugo2016_S2B.csv` | mmc2 S2B | 72×7 | **Differentially enriched gene sets** (the IPRES components): p/FDR/stat/diffAvg. |
| `hugo2016_S2C.csv` | mmc2 S2C | 19×3 | **IPRES gene-set definitions**: geneset name, detail, member gene listing. |

## Use in this project

- **S2B/S2C → IPRES signature**: the published innate-anti-PD-1-resistance
  signature. We score it as a benchmark our latent RNA-state axis must beat or
  explain (EMT, angiogenesis, wound-healing, hypoxia programs).
- **S1D MAF + S1B → driver mutations & true TMB**: the confounder-to-beat and
  the driver-origin material (BRAF/NRAS/BRCA2), finer than cBioPortal's
  harmonized 27-sample release.
- **S1B Purity/Ploidy → composition/instability covariates.**
- **S1C/S1B neoepitopes → antigenicity endpoint.**

## Verification (against the paper)

- WES median coverage 140× — `S1B.AvgCov` median = 140.3 ✅
- median 489 non-synonymous mutations — `S1B.TotalNonSyn` median = 489 (n=38, range 73–3985) ✅
- BRCA2 responder-enrichment — cross-tab of `S1B.BRCA2mutant?` × Response:
  6 R / 1 NR among BRCA2-mutant vs 15 R / 16 NR among WT (directional
  enrichment confirmed; formal significance not tested here) ✅
- IPRES = MAPKi-induced EMT/angiogenesis/wound-response gene sets — present in S2B/S2C ✅
- Response classes — S1A irRECIST CR 7 / PR 14 / PD 17 (38 patients with a
  non-null irRECIST call in the 42-row S1A table; **not** a WES-defined subset) ✅

Header rows in the raw workbooks start at row 3 (0-indexed 2); title and blank
rows above were dropped during extraction. Reproduce via
`clean/_extract_manifest.json`.
