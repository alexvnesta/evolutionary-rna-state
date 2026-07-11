# Verified Data Inventory — 2026-07-11 (single source of truth)

Computed fresh from disk + `results/predictor/phase2_covariates_n106.parquet` (the canonical label/cohort table).
Reconciles every count that has appeared in discussion. Machine-readable: `results/DATA_INVENTORY.csv`.

## Why counts differed before
Numbers vary **by layer**, not by error. Three distinct things were being counted interchangeably:
1. **Physical aligned BAMs** (40) — the alignment substrate.
2. **Per-layer caller outputs** (junctions 59, TE 41, IR 40, AEI 43) — each layer has different coverage
   because outputs accumulated across several runs + an older 16-sample `editing_bams` set.
3. **The canonical labeled cohort** (106) — the target denominator from the covariate table.

## LAYER 0 — physical aligned BAMs on disk: **40**
| set | count | note |
|---|---|---|
| hisat2 markdup (`results/rnaseq_cohort/hisat2/`) | 40 | **30 Gide + 10 Riaz** — the current caller substrate |
| STAR (`results/rnasplice_cohort/star_salmon/`) | 40 | SAME 40 samples, different aligner (rnasplice) |
| editing_bams (`results/editing_bams/`, ERR-named) | 16 | older separate set, partly overlapping |

The 40 hisat2 BAMs are **30 Gide + 10 Riaz + 0 Hugo**. (Hugo junction/TE/IR outputs come from the older
editing_bams + a sibling's separate editing run, not these 40.)

## LAYER 1 — feature coverage across the 106 labeled samples
| cohort | n labeled | immune floor | junctions (splice) | TE | IR | AEI (editing) |
|---|---:|---:|---:|---:|---:|---:|
| gide2019 | 69 | 57 | 32 | 14 | 13 | 14 |
| hugo2016 | 27 | 27 | 22 | 22 | 22 | 24 |
| riaz2017 | 10 | 10 | 5 | 5 | 5 | 5 |
| **TOTAL** | **106** | **94** | **59** | **41** | **40** | **43** |

## What this means for the analysis
- **Immune floor (expression)** is the most complete layer (94/106) — it needs only quantification, not the full caller stack.
- **Splicing/junctions** (59) is the layer the Evo2 encoder test uses — best in Gide (32) and Hugo (22).
- **The two VALID within-cohort test frames** (need floor + a non-ref layer + a working positive control):
  - **Gide n=32** (junctions) — the best-powered frame; Evo2 test done, p=0.83, null.
  - **Hugo n=22** (junctions) — replication; Evo2 test done, p=0.259, null.
  - Riaz n=5 — too small for a within-cohort test.
- **Cross-cohort LOCO is INVALID** for this 3-cohort set (the immune floor itself fails to transfer, Phase 2), so
  every test is within-cohort.

## The gap to a fully-powered result
To deepen **Gide 32→69** (double the best frame) requires aligning **37 more Gide samples**. There are **NO source
FASTQs on local disk** (deleted post-alignment; only HLA fragments remain), so this means **re-downloading ~300 GB
from ENA**. Local disk (399 GB free) + saturated local CPU (a sibling runs AEI+rnasplice at load ~18) make local
alignment impractical now — hence alignment is being moved to **Modal** (off-local disk + CPU).
