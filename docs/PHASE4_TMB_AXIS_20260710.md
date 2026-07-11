# Phase 4 — is a burden (TMB) axis more cross-cohort-portable than the inflammation floor?

**Date:** 2026-07-10 | **Status:** COMPLETE (exploratory, underpowered) | **Compute:** local, no GPU.

## Motivation
Phase 3 showed the inflammation floor's failure to transfer is because its response DIRECTION flips across
cohorts (Hugo's response is not inflammation-driven; only low TGFβ tracks it there). Question: is there ANY
axis that points the same way across cohorts? TMB is the obvious candidate. Phase 2's audit correction
established TMB IS on disk for Hugo (25/27) and Riaz (10/10), just under patient-level ids.

## TMB join (recovered the crosswalk the earlier phase missed)
`run_catalog.csv` maps `run_accession → patient_id`; TMB is keyed by patient (Hugo `Pt#`, Riaz `Pt#_pre`).
Joined: **Hugo 25/27, Riaz 10/10, Gide 0/69** (no Gide TMB exists on disk — genuinely absent). Merged into
`results/predictor/phase2_covariates_n106.parquet`.

## Result — TMB is directionally consistent, unlike the floor (but underpowered)
| | Hugo (n=25) | Riaz (n=10) |
|---|---|---|
| TMB→response correlation | r=+0.25 (p=0.227) | r=+0.28 (p=0.438) |

Both **positive** — TMB points the same way in both cohorts. Contrast with the inflammation floor, whose
GEP correlation was +0.59 (Gide) / −0.05 (Hugo) / +0.24 (Riaz) — a sign flip.

Cross-cohort transfer (single TMB feature): **Hugo→Riaz 0.520, Riaz→Hugo 0.635**. Weak, but NOT the
below-chance collapse the full floor showed (0.200–0.520 pairwise). Within-Hugo, TMB does not add to the
floor (floor 0.689 → floor+TMB 0.679; TMB-alone 0.564).

## Honest caveats
- **Underpowered:** n=35 total with TMB, neither within-cohort correlation is significant (p=0.23, 0.44).
  The "directional consistency" is a sign agreement on noisy small samples, not a demonstrated effect.
- **No Gide TMB:** the largest cohort has none, so a 3-cohort TMB LOCO is impossible; this is a 2-cohort
  (Hugo/Riaz) observation only.
- TMB does not RESCUE cross-cohort prediction to a useful level (0.52–0.64 is still weak).

## What this contributes to the program
1. It identifies the SHAPE of a portable signal: a burden/genomic-instability axis is directionally stable
   across cohorts where inflammation is not. This is a concrete hypothesis for what a cross-cohort-valid
   feature must look like — and it happens to be exactly the axis the sequence-model (Evo 2) branch is meant
   to probe (aberrant-sequence load ≈ a burden-like quantity), reinforcing that the encoder's target should
   be burden/aberrancy, not inflammation re-derivation.
2. It does NOT change the gate: cross-cohort inference remains underpowered and Gide-less for TMB. Within-
   Gide stays the only frame with a strong positive control.
3. Actionable: acquiring Gide TMB (or computing a burden proxy from the RNA, e.g. expressed-mutation load)
   is a higher-value, lower-cost next step than the GPU encoder pass, because it directly tests whether the
   directionally-consistent axis holds at 3-cohort scale.
