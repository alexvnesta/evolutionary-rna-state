# Phase 3 — transferability investigation: is the LOCO collapse a regimen confound?

**Date:** 2026-07-10 | **Status:** COMPLETE | **Compute:** local, no GPU.
**Question:** Phase 2 found the immune floor fails cross-cohort (LOCO 0.507 vs within-cohort 0.68–0.81),
with the floor→response direction flipping across cohorts. Is that driven by REGIMEN heterogeneity
(Gide mono+combo; Riaz anti-CTLA4-progressed) — i.e. removable by matching treatment arm?

## Regimen annotation (106/106)
| cohort | anti-PD1 mono | anti-PD1+CTLA4 combo | anti-PD1 (Riaz, mixed prior-ipi) |
|---|---|---|---|
| gide2019 | 38 | 31 | – |
| hugo2016 | 27 | – | – |
| riaz2017 | – | – | 10 |
(Gide arm from `data/registry/gide_id_crosswalk.csv`; Hugo/Riaz cohort-level.)

## Answer: NO — regimen does NOT explain the failure
Restricting to the SAME regimen (anti-PD1 monotherapy: Gide-mono n=32 vs Hugo n=27):
- Per-cohort floor→response correlation STILL flips: Gide-mono GEP **+0.68**, Hugo GEP **−0.05**.
- Matched-regimen cross-cohort transfer stays at chance: **Gide-mono→Hugo 0.556, Hugo→Gide-mono 0.576.**

So the direction-flip survives regimen matching. It is a cohort/biology/label property, not a treatment-arm
confound.

## Mechanism: Hugo's ICB biology genuinely differs from Gide's
Within Hugo (n=27), response separation per floor feature (Mann-Whitney R vs NR):
| feature | p | responder vs non-responder median |
|---|---|---|
| gep_tcell_inflamed | 0.75 | 0.06 vs 0.37 (WRONG direction) |
| ifng_score | 0.79 | 0.07 vs 0.32 (WRONG direction) |
| teff | 0.79 | −0.05 vs 0.44 (WRONG direction) |
| **tgfb** | **0.02** | −0.53 vs 0.72 (responders LOW TGFβ) |
| teff_tgfb_balance | 0.08 | 0.44 vs −0.57 |

In Gide the T-cell-inflamed GEP is the dominant positive predictor; **in Hugo it is non-significant and
points the wrong way**, and only low TGFβ tracks response. The cohorts encode different response biology
(consistent with Hugo's known smaller n, different sequencing/label pipeline, and possible occult
composition differences). A single cross-cohort decision boundary cannot fit both.

## Firm consequences (updates the plan)
1. **LOCO is not a valid frame for this 3-cohort set** — confirmed a second way (regimen-matched transfer
   also fails). Cross-cohort claims are off the table until more/harmonized cohorts exist.
2. **Within-cohort (Gide) is the only frame with a working positive control** — any encoder test must live
   here, reported descriptively, with cross-cohort generalization stated as untestable. Gide is 69 samples
   total; the floor positive control runs on the 57 with complete floor features + labels (n=57 is the
   analysis subset, NOT the cohort size).
3. **TMB is worth adding for Hugo/Riaz** (crosswalk-joinable; see Phase 2 doc) since inflammation is NOT
   Hugo's predictor — a TMB/burden axis may be what transfers there. This is a concrete no-GPU follow-up.
4. The encoder GPU pass remains correctly deferred: the frame it would be judged in does not support
   cross-cohort inference, and within-Gide is the tractable target once the sequence pipeline is built.
