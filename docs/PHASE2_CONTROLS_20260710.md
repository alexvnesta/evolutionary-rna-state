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
- **TMB: UNAVAILABLE for this cohort.** The on-disk `tmb_standardized.parquet` is keyed to dfci2019/liu2019
  WES cohorts (generic `Sample*` ids) with ZERO overlap with our RNA accessions. There is no crosswalk. TMB
  as a competing covariate cannot be assembled without new data. Documented, not silently skipped.
- **Purity:** available for the 21–25 non-ref samples (prior InstaPrism frame).

## Recommendation
Do NOT proceed to a cross-cohort encoder scoring pass until the floor-transfer failure is addressed. The
highest-value next work is batch harmonization + a re-test of floor LOCO transfer — a no-GPU experiment that
gates everything downstream. Within-cohort (Gide, n=57) analyses remain valid and are where any encoder
signal should first be sought.
