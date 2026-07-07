# Riaz et al. 2017 — supplemental data

**Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab.**
Riaz N, et al. *Cell* 171(4):934–949.e16, 2017. PMID 29033130.
Nivolumab (Ipilimumab-naive + Ipilimumab-progressed); longitudinal
WES + RNA-seq + TCR-seq, pre- and on-treatment.

Source files (publisher supplemental): `mmc1.xlsx`–`mmc6.xlsx`.
MD5s in `raw/checksums.json`.

## Cleaned tables (`clean/`)

| CSV | Source | Rows×Cols | Contents |
|-----|--------|-----------|----------|
| `riaz2017_S1.csv` | mmc1 | 68×40 | Exome sequencing QC per sample: total reads, mean target coverage, PCT_TARGET_BASES at 20–100×. |
| `riaz2017_S2.csv` | mmc2 | 73×12 | **Clinical + summary genomics per patient**: Cohort, Response, survival, subtype, mutational subtype, M-stage, **Mutation Load, Neo-antigen Load, Neo-peptide Load, Cytolytic Score**. |
| `riaz2017_S3.csv` | mmc3 | **28589×11** | **Somatic mutations (MAF)**: Patient, Hugo Symbol, Chr, Start/End, Variant Classification, HGVS c/p, **Tcov, Tac, Taf** (tumor coverage / alt count / allele frequency). 68 patients. |
| `riaz2017_S4.csv` | mmc4 | 80×7 | **Genomic data-availability matrix**: which of {Exome, TCR-seq, RNA-seq} exist pre/on-treatment per patient. |
| `riaz2017_S5A.csv` | mmc5 | 41477×8 | Class I neoantigens (set A): WT/MT peptide, WT/MT allele, WT/MT binding score. |
| `riaz2017_S5B.csv` | mmc5 | 4492×8 | Class I neoantigens (set B, filtered). |
| `riaz2017_S6A–D.csv` | mmc6 | 189 / 1384 / 475 / 2670 ×6 | Differential expression contrasts: symbol, baseMean, log2FC, lfcSE, pvalue, padj. |

## Use in this project

- **S3 MAF with Taf → clonality/CCF**: allele frequencies + the paper's
  PyClone/FACETS framework let us anchor clonal-vs-subclonal structure — the
  clonal-evolution arm. (Per-SNV CCF can be recomputed from Taf × purity.)
- **S2 → true TMB, neoantigen load, cytolytic score**: confounders-to-beat and
  the antigenicity endpoint, per patient.
- **S4 → longitudinal design**: identifies patients with paired pre/on samples
  across modalities — the "tumor evolves along a trajectory" test set.
- **S5A/B → antigenicity** (Class I neoantigen catalog).
- **S6A–D → response-associated DE programs** for benchmarking.

## Verification (against the paper)

- 73 clinical patients (S2) — matches iAtlas harmonized cohort ✅
- Mutation-load median 182 (paper reports 183; range 1–7360 exact) ✅
- 68 patients with somatic mutation calls (S3) ✅
- Longitudinal multimodal design (WES+RNA+TCR, pre/on) present in S4 ✅

Header rows differ per sheet (title + notes above the header); the extractor
auto-detects the header row. Newlines inside column names were flattened.
Reproduce via `clean/_extract_manifest.json`.
