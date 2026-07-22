# evolutionary-rna-state

**Reconstructing a latent evolutionary RNA-state of tumors from bulk RNA-seq to organize immune checkpoint blockade (ICB) response.**

![The evolutionary RNA-state hypothesis: early driver mutations → genomic instability + epigenetic remodeling → RNA-processing dysregulation → latent RNA-state S → tumor antigenicity and ICB response](results/cartoon_thesis.png)

*Early driver mutations set a tumor's trajectory; downstream RNA-processing phenotypes are manifestations of one latent state **S** that shapes antigenicity and immunotherapy response. Illustrations: [NIH BioArt](https://bioart.niaid.nih.gov) (NIAID Visual & Medical Arts).*

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

## Status: ARCHIVED (2026-07-22) — rigorous negative with mechanism

**Start here if reopening:** [`docs/ARCHIVE_SUMMARY_20260722.md`](docs/ARCHIVE_SUMMARY_20260722.md)
— one-page re-entry point (verdict, what's dead, the one live path, where everything lives).

**Final verdict.** The strong form of the hypothesis is **falsified on bulk RNA-seq**: no
tumor-intrinsic RNA representation carries ICB-response information beyond immune composition,
and the immune signal itself does not transport across cohorts. The response-tracking signal in
bulk tumor RNA *is* immune composition. This is a cautionary negative **with a mechanism**: the
tumor-intrinsic aberrancy signal is a downstream readout of the inflamed state whose coupling to
inflammation **sign-flips across cohorts**, which is why cross-cohort transfer (LOCO) fails.

**Diligence is complete on bulk.** Tested and null: WES antigen quantity (416 samples, perm
p=0.78); every non-reference RNA burden layer; four frozen sequence encoders (EVA, HyenaDNA,
Evo2, Caduceus — 4/4); trained MLP and deep VAE (scVI); a structure-aware dsRNA/viral-mimicry
feature (novel: decoupled from interferon, yet still null); the literal latent-state-vs-clonality
claim (first powered test, n≈60); and the domain-invariance method class (IRM/GroupDRO) built
specifically for the sign-flip. The remaining open path is **substrate, not method**: a matched
single-cell malignant compartment, not more architecture search.

**Dated deliverables (2026-07-22):**

| Deliverable | Path |
|---|---|
| **Archive summary (re-entry point)** | [`docs/ARCHIVE_SUMMARY_20260722.md`](docs/ARCHIVE_SUMMARY_20260722.md) |
| Deep-learning diligence ledger | [`docs/DEEP_LEARNING_DILIGENCE_20260722.md`](docs/DEEP_LEARNING_DILIGENCE_20260722.md) |
| dsRNA viral-mimicry + latent-state results | [`docs/RESULTS_VIRALMIMICRY_LATENTSTATE_20260722.md`](docs/RESULTS_VIRALMIMICRY_LATENTSTATE_20260722.md) |
| Gemini-exchange assessment | [`docs/GEMINI_ASSESSMENT_20260722.md`](docs/GEMINI_ASSESSMENT_20260722.md) |
| Pre-registered fine-tune design (next step) | [`docs/EXPERIMENT_SCOPE_FINETUNE_RAWREADS_20260722.md`](docs/EXPERIMENT_SCOPE_FINETUNE_RAWREADS_20260722.md) |

<details><summary>Earlier submission deliverables (pre-archival, superseded headline)</summary>

An earlier de-novo antigen-presentation axis separated responders *within* a cohort (LOO AUROC
0.87) but was infiltration-driven (ρ=0.77) and did not transfer (Riaz 0.36, Hugo 0.58) — the
first signature of the non-transportability now explained mechanistically above.

| Deliverable | Path |
|---|---|
| Manuscript (Markdown) | [`docs/WRITEUP.md`](docs/WRITEUP.md) |
| **Manuscript (HTML, linked citations + cartoons)** | [`docs/manuscript.html`](docs/manuscript.html) |
| Mechanistic cartoons (NIH BioArt illustrations) | [`results/cartoon_thesis.png`](results/cartoon_thesis.png), `cartoon_wes_blind.png`, `cartoon_pipeline.png`, `cartoon_confound.png` |
| BioArt icon set + attribution | [`results/bioart_icons/`](results/bioart_icons/) (NIAID Visual & Medical Arts) |
| References (verified DOIs) | [`results/citations.json`](results/citations.json) |
| 100–200 word summary | [`docs/SUMMARY_100-200w.md`](docs/SUMMARY_100-200w.md) |
| Figure deck (8 figs) + supplement | [`results/figure_deck.pdf`](results/figure_deck.pdf), `results/figure_deck_supplement.pdf` |
| Reproducible demo (< 30 s) | [`notebooks/demo_reproduce_headline.ipynb`](notebooks/demo_reproduce_headline.ipynb) |
| How to reproduce everything | [`REPRODUCE.md`](REPRODUCE.md) |
| Demo video script + slideshow | [`docs/VIDEO_SCRIPT.md`](docs/VIDEO_SCRIPT.md), `results/demo_slideshow.gif` |
| Daily-refresh skill | `evolutionary-rna-state-refresh` (published; trigger with "run the refresh") |
| BioArt figure skill | `niaid-bioart` (published; search + download NIH BioArt illustrations for figures) |

</details>

## Provenance & compliance

This repository was initialized fresh during the hackathon; all code is
authored in-session. Public-data acquisition is permitted and the committed
`data/` payload is metadata for already-public datasets. See `COMPLIANCE.md`.

## License

MIT — see `LICENSE`.
