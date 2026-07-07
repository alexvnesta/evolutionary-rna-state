# Gide et al. 2019 — supplemental data

**Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and
Anti-PD-1/Anti-CTLA-4 Combined Therapy.** Gide TN, et al. *Cancer Cell*
35(2):238–255.e6, 2019. PMID 30753825. RNA-seq: ENA PRJEB23709.
Anti-PD-1 monotherapy and combined anti-PD-1/anti-CTLA-4.

Source files (publisher supplemental): `mmc2.xlsx`–`mmc9.xlsx` (Tables S1–S8),
`mmc1.pdf` and `mmc10.pdf` (supplemental figures/legends).
MD5s in `raw/checksums.json`. **This cohort has no WES** — supplements are
clinical + transcriptomic only, consistent with the cBioPortal record.

## Cleaned tables (`clean/`)

| CSV | Source | Rows×Cols | Contents |
|-----|--------|-----------|----------|
| `gide2019_S1_PD1_clinical.csv` | mmc2 S1 | 54×15 | **Anti-PD-1 monotherapy patients**: age, sex, treatment, Best RECIST, PFS, OS, PRE/EDT biopsy sites, timepoints. |
| `gide2019_S2_comboclinical.csv` | mmc3 S2 | 51×15 | **Combined anti-PD-1/anti-CTLA-4 patients**: same fields. |
| `gide2019_S3_PD1_DEgenes.csv` | mmc4 S3 | 310×14 | DE genes R vs NR (mono, adj p<0.05): gene, locus, adj p, cluster, description, per-group expression. |
| `gide2019_S4_PD1_KEGG.csv` | mmc5 S4 | 12×7 | KEGG gene sets enriched R vs NR (mono): SIZE, ES, NES, NOM p, FDR q, DE genes. |
| `gide2019_S5_combo_DEgenes.csv` | mmc6 S5 | 328×13 | DE genes R vs NR (combo). |
| `gide2019_S6_combo_KEGG.csv` | mmc7 S6 | 11×7 | KEGG gene sets enriched R vs NR (combo). |
| `gide2019_S7_Venn_PRE.csv` | mmc8 S7 | 123×3 | DE-gene overlap pre-treatment (R vs NR). |
| `gide2019_S8_Venn_EDT.csv` | mmc9 S8 | 133×3 | DE-gene overlap early-during-treatment. |

`mmc1.pdf` / `mmc10.pdf`: supplemental figures (mass-cytometry gating,
EOMES+CD69+CD45RO+ effector-memory T-cell analyses, additional expression
panels). Kept in `raw/`; parse specific pages on demand.

## Use in this project

- **S1/S2 → clinical backbone + arm stratification**: the mono-vs-combo split
  is a within-cohort confounder we must stratify or model; these tables carry
  the per-patient arm, RECIST, PFS/OS, and biopsy timepoints.
- **S3–S6 → response-associated DE programs & KEGG pathways** for both arms:
  benchmarks the latent RNA-state axis competes against.
- **S7/S8 → PRE vs EDT program overlap**: a longitudinal-consistency check.
- Gide is the **RNA development lead** (deepest/longest reads, balanced
  response); its value is the RNA substrate + immune-composition benchmarks,
  not genomics.

## Verification (against the paper)

- Two treatment arms (anti-PD-1 mono; anti-PD-1/anti-CTLA-4 combo) — S1/S2 ✅
- Clinical tables are **patient-level**: 54 mono + 51 combo = 105 patients
  (verified: mmc2 S1 contains exactly 54 non-null patient rows numbered 1–54, so
  extraction was not clipped). The paper abstract's "n=63 / n=57" are larger than
  these patient counts; the most likely explanation is that those are
  biopsy/sample-level counts (the study spans multiple biopsy timepoints), but
  **this has not been confirmed against the Gide 2019 main text or a biopsy
  manifest in-session** — treat the patient-vs-sample reconciliation as a
  hypothesis, not a settled fact.
- RNA-seq accession **PRJEB23709** — verified via ENA: study title "Biomarkers
  of response and resistance to checkpoint blockade immunotherapy in metastatic
  melanoma", RNA-Seq, Homo sapiens.
- PMID **30753825** — verified via EuropePMC lookup on DOI 10.1016/j.ccell.2019.01.003.
- RECIST/PFS/OS present per patient ✅

Header row auto-detected per sheet; column-name newlines flattened. Reproduce
via `clean/_extract_manifest.json`.
