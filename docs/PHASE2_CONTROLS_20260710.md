# Phase 2 — covariates + biological positive controls (the gate)

**Date:** 2026-07-10 | **Status:** COMPLETE, with one decisive negative finding | **Compute:** local, no GPU.

## Purpose
Before spending GPU on an encoder, the protocol (§6) requires biological positive controls to PASS — chiefly
that the immune floor recovers ICB response under the SAME leave-one-cohort-out (LOCO) frame the encoder
would be judged by. If the positive control cannot transfer cross-cohort, no downstream feature can be
credited with cross-cohort signal either.

## Headline finding — the immune floor does NOT transfer cross-cohort at cohort scale
| Frame | AUROC | n |
|---|---|---|
| Within-Gide 5-fold | 0.807 | 57 |
| Within-Hugo 5-fold | 0.689 | 27 |
| Within-Riaz 5-fold | 0.680 | 10 |
| **LOCO 3-cohort (full 5-feature floor)** | **0.507** | 94 |
| GEP-alone LOCO | 0.621 | 94 |
| teff/tgfb-balance-alone LOCO | 0.666 | 94 |

Every pairwise cross-cohort transfer is at or below chance: gide→hugo 0.500, gide→riaz 0.520, hugo→gide
0.515, hugo→riaz 0.200, riaz→gide 0.427, riaz→hugo 0.233.

## This is real, not an artifact — three checks
1. **Not label imbalance:** response rates are balanced across cohorts (Gide 0.579, Hugo 0.556, Riaz 0.500).
2. **Signal IS present within each cohort** (0.68–0.81), so the features are informative — they just don't
   transfer.
3. **Batch structure dominates:** cohort identity is predictable from the floor features at acc=0.606 (chance
   0.333) — the feature space carries more cohort/batch signal than transferable biology at this n and
   harmonization level.

## Consequence for the encoder phase (decisive)
The positive control **FAILS the LOCO gate at n≈94**. Per the protocol's own logic, this means:
- A 3-cohort LOCO encoder test is **not currently interpretable** — a null encoder result would be
  unattributable (could be "no signal" or "same cross-cohort collapse the floor shows"), and a positive
  result would be suspect.
- The binding constraint is **NOT the model** (Evo 2 vs alternatives) — it is **cross-cohort
  transferability / batch harmonization at achievable n**. Adding a sequence model cannot fix a floor that
  itself does not transfer.
- This SHARPENS the descriptive-regime decision (§0a): not only is significance unreachable with 3 cohorts,
  the cross-cohort *point estimate* itself is dominated by batch effects. Cross-cohort claims require either
  (a) more cohorts, or (b) explicit batch harmonization (e.g. ComBat/rank-normalization within the floor
  space) validated to restore floor transfer BEFORE any encoder is added.

## Covariate availability (honest inventory)
- **Immune floor:** 94/106 samples. Available.
- **HLA-I het/LOH:** het fraction for the 16 non-ref-matrix samples (mean 0.885). Available (this session).
- **TMB: present for 2/3 cohorts but not yet joined (ID-crosswalk gap, NOT missing data).** The on-disk
  `tmb_standardized.parquet` (336 rows) covers dfci2019 (144) and liu2019 (122) — unrelated WES cohorts — but
  ALSO carries **hugo2016 (26, real tmb_rate values, keyed `Pt1`-style) and riaz2017 (44, keyed `Pt10_pre`)**.
  It has **NO gide2019 rows** (0), i.e. no TMB for our largest cohort. The Hugo/Riaz TMB is keyed by
  patient-level ids, not `run_accession`, so a direct accession join gives zero overlap — but those ids are
  joinable via the existing crosswalks (`data/registry/`). So TMB is assemblable for Hugo+Riaz with a
  crosswalk step; it is genuinely absent only for Gide. This does not change the Phase-2 gate conclusion
  (the floor-transfer failure is independent of TMB), but TMB is NOT "unavailable without new data" — that
  earlier phrasing was wrong.
- **Purity:** available for the 21–25 non-ref samples (prior InstaPrism frame).

## Recommendation
Do NOT proceed to a cross-cohort encoder scoring pass until the floor-transfer failure is addressed. The
highest-value next work is batch harmonization + a re-test of floor LOCO transfer — a no-GPU experiment that
gates everything downstream. Within-cohort Gide analyses (n=57 floor-complete of 69) remain valid and are where any encoder
signal should first be sought.

## ROOT CAUSE (Phase 2b) — the association DIRECTION flips across cohorts
Harmonization (within-cohort rank-normalization, within-cohort z-scoring) does NOT restore transfer
(0.507 → 0.531 / 0.519). Reason, from per-cohort point-biserial correlation of each floor feature with
response:

| feature | Gide (n=57) | Hugo (n=27) | Riaz (n=10) |
|---|---|---|---|
| gep_tcell_inflamed | +0.59 | −0.05 | +0.24 |
| ifng_score | +0.60 | −0.05 | +0.13 |
| teff | +0.58 | −0.08 | +0.15 |
| tgfb | −0.02 | **−0.50** | **+0.61** |
| teff_tgfb_balance | +0.43 | +0.31 | −0.44 |

The canonical inflammation features (GEP/IFNG/teff) are strong positive predictors in Gide, ~null in Hugo,
weak in Riaz; `tgfb` points in OPPOSITE directions in Hugo vs Riaz. A classifier trained on one cohort's
direction is miscalibrated or reversed on another. **This is not removable by feature-space batch correction**
— it reflects genuine cohort heterogeneity (regimen mix: Gide mono+combo, Riaz anti-CTLA4-progressed;
biopsy/label differences; small-n noise in Hugo/Riaz). It is the mechanistic reason LOCO collapses.

**Firm consequence:** cross-cohort (LOCO) is not a valid evaluation frame for THIS 3-cohort set at THIS n,
for ANY feature block including a sequence encoder. Within-cohort Gide (n=57 floor-complete of 69 total) is the only frame where a
positive control holds. The encoder hypothesis, if pursued now, must be tested WITHIN Gide with the honest
caveat that cross-cohort generalization is untestable until the cohort-heterogeneity/regimen confounds are
resolved (more cohorts, harmonized regimens, or per-regimen stratification).
