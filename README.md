# evolutionary-rna-state

**Reconstructing a latent evolutionary RNA-state of tumors from bulk RNA-seq to organize immune checkpoint blockade (ICB) response.**

## Thesis

Early driver mutations set a tumor's evolutionary trajectory. As the tumor
evolves, genomic instability, epigenetic remodeling, RNA-processing
dysregulation, and immune selection combine to produce coordinated
transcriptomic abnormalities — alternative splicing, intron retention, RNA
editing, transposable-element (TE) activation, fusion transcripts, and
cryptic/non-canonical ORFs.

**Core reframing:** these are not independent biomarkers but downstream
manifestations of a single latent *evolutionary RNA-state* **S**. That state —
not any one biomarker — ultimately shapes tumor antigenicity and response to
ICB. Clinical response is one noisy *observable* of where a tumor sits on its
trajectory, not the latent variable itself.

## Falsifiable claims

1. **Co-variation (internal).** The RNA phenotypes above share variance — a
   low-rank structure exists — rather than behaving independently. Tested
   *without reference to response labels.*
2. **Organization (external).** A low-dimensional representation built to
   capture that shared RNA-state variance also stratifies ICB response, and
   does so beyond the field-standard confounders (TMB / expressed-neoantigen
   load, tumor purity, and immune/stromal composition).

## Why bulk RNA-seq, why raw reads

Annotation-based pipelines discard non-reference signal — which is exactly
where evolutionary RNA-state fingerprints live. The design pairs interpretable
expression/signature features with a raw-read encoder branch so both
reference and non-reference signal can contribute to the sample representation.

## Repository layout

| Path | Contents |
|------|----------|
| `data/` | **Metadata only** — manifests, download receipts, run catalog, clinical/response labels. No raw reads (see `data/README.md`). |
| `src/` | Library code (authored in-session). |
| `analysis/` | Analysis scripts / pipeline stages. |
| `notebooks/` | Exploratory notebooks. |
| `results/` | Generated tables and figures (large outputs git-ignored). |
| `docs/` | Data inventory, roadmap, methods notes. |

## Data at a glance

Development/validation centers on pretreatment melanoma ICB RNA-seq with public
raw reads, with additional public IO cohorts inventoried for extension:

- **Gide 2019** (melanoma, anti-PD-1 ± anti-CTLA-4) — ENA `PRJEB23709`
- **Riaz 2017** (melanoma, nivolumab) — SRA/ENA `PRJNA356761` / GEO `GSE91061`
- Additional inventoried cohorts: Zhao/Cloughesy 2019 (GBM), Kim 2018
  (gastric), IMvigor210 (urothelial, processed), TISCH2 (scRNA validation).

See `docs/DATA_INVENTORY.md` and `data/dataset_summary.csv` for the full table,
access levels, and accessions.

## Provenance & compliance

This repository was initialized fresh during the hackathon; all code is
authored in-session. Public-data acquisition is permitted and the committed
`data/` payload is metadata for already-public datasets. See `COMPLIANCE.md`.

## License

MIT — see `LICENSE`.
