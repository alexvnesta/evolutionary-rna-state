# Reproducing the results

This project reconstructs a candidate latent "evolutionary RNA-state" of tumors
from bulk RNA-seq and tests whether it organizes immune-checkpoint-blockade
(ICB) response beyond established biomarkers. All data are public; all code was
authored during the event.

## Quick start (< 30 s, no raw-read download)

The headline result reproduces from committed per-sample gene-TPM matrices:

```bash
pip install -r requirements.txt
jupyter notebook notebooks/demo_reproduce_headline.ipynb   # or: run the 4 cells
```

The notebook loads the de-novo salmon-quantified gene TPM for 52 pre-treatment
melanoma transcriptomes (Gide n=30, Riaz n=10, Hugo n=12), builds the
antigen-presentation axis, and reproduces:

| Result | Value |
|---|---|
| Within-Gide fold-contained LOO AUROC (B2M/HLA-A/TAP1) | **0.888** |
| Antigen axis ↔ immune-infiltration axis (Spearman ρ) | **0.77** |
| Gide→Riaz held-out AUROC | **0.36** |
| Gide→Hugo held-out AUROC | **0.58** |

**Read:** a de-novo antigen axis that separates responders within one cohort is
largely infiltration-driven and does **not** transfer to either independent
held-out cohort — a rigorous, cautionary negative.

## Full analysis

| Stage | Code | Output |
|---|---|---|
| Data freeze / harmonization (5 cohorts) | `src/data.py` | `results/analysis_frame.parquet` |
| Internal co-variation null (WES proxies) | `analysis/01_covariation.py` | `results/covariation_*`, `fig_covariation.png` |
| Latent-state + response modeling machinery | `src/model.py` | `results/model_validation.json` |
| Raw-read de-novo quantification (salmon) | `analysis/pilot/run_salmon_pilot.sh` | `results/*_salmon/` (FASTQ streamed + deleted) |
| Per-sample gene TPM assembly | `analysis/pilot/analyze_pilot.py` | `results/features/quant_gene_tpm.parquet` |
| Deepened de-novo analysis (n=40) | `analysis/pilot/deepen_analysis.py` | `results/deepen_*`, `fig_deepen.png` |
| Three-cohort held-out transfer | (this session) | `results/heldout_transfer_3cohort.json`, `fig_transfer_3cohort.png` |
| Rigor + robustness | (this session + sibling harness) | `results/denovo_robustness.json`, `fig_robustness.png` |

## Raw-read pipeline (optional, hours)

The raw-read arm streams the first 3M read-pairs per sample from ENA, quantifies
with salmon against a GENCODE v44 transcriptome index, verifies gzip integrity,
and deletes the FASTQ (stream-align-delete; ~700 GB free < 1.3 TB full reads):

```bash
# build index once
salmon index -t refs/gencode.v44.transcripts.fa.gz -i refs/gencode_v44_index --gencode -k 31 -p 8
# quantify a manifest (cohort,run_accession,...,fastq_ftp,fastq_md5)
NREADS=3000000 THREADS=8 GENEMAP=refs/tx2gene.tsv \
  bash analysis/pilot/run_salmon_pilot.sh <manifest.csv> refs/gencode_v44_index <out_dir>
```

## Deliverables

- **Manuscript:** `docs/WRITEUP.md`
- **Figure deck (8 figures):** `results/figure_deck.pdf` (+ `results/figure_deck_supplement.pdf`)
- **100–200 word summary:** `docs/SUMMARY_100-200w.md`
- **Reproducible demo:** `notebooks/demo_reproduce_headline.ipynb`

## Compliance

Repository initialized at event kickoff; all code authored in-session. Only
public data used (cBioPortal/iAtlas metadata + ENA FASTQ). Raw reads are
git-ignored and never redistributed. See `COMPLIANCE.md`.
