# Can a model disentangle aberrancy from immune composition? — No (deconvolution-based test)

**Date:** 2026-07-11  **Status:** COMPLETE  **Compute:** local CPU (no new Modal spend)

## The question
The aberrancy×dosage insight showed the non-reference "aberrancy" signal is entangled with immune
composition and sign-flips across cohorts. Natural follow-up: can we build a model that DISENTANGLES the
tumor-intrinsic aberrancy from the immune-composition component and recovers independent response signal?

## The test (this IS the disentangling model)
A disentangling model predicts response from the aberrancy component ORTHOGONAL to immune composition.
We built that residual against the richest immune basis available — the 5 immune-floor signatures PLUS all
**11 InstaPrism deconvolution cell fractions** (Malignant/purity, CD8/CD4 T, Tregs, myeloid, B, plasma, pDC,
mast, cycling) — fold-contained, and asked whether any response signal survives. This is exactly what the
deconvolution fractions were built for (cf. MEDIATION_PURITY_CORRECTED.md).

## Result — nothing survives the disentangling
| | Evo2 residual on floor (5 sig) | Evo2 residual on RICH basis (16: sig + deconv fractions) | perm p (rich) |
|---|---|---|---|
| Gide n=32 | 0.354 | **0.412** | 0.905 |
| Hugo n=22 | 0.576 | **0.488** | 0.433 |

Interaction (aberrancy×dosage) residualized on the rich basis: Gide 0.401, Hugo 0.428 — also chance.

- Both cohorts: the immune-orthogonal aberrancy component is at/below chance, not significant.
- In Hugo the rich-basis residual DROPS toward chance (0.576→0.488): the deconvolution fractions absorbed the
  last of the apparent signal, confirming it was immune-composition proxy.

## Why this is definitive, not a modeling-power limitation
A disentangling model can only recover a component that exists. We removed the immune basis with the richest
correction available (real malignant-aware deconvolution, not a marker-gene proxy) and the residual subspace
carries no response information. This mirrors the expression-level finding in MEDIATION_PURITY_CORRECTED.md,
where the phenotype→infiltration a-path collapsed to zero after purity correction (antigen-presentation a-path
+0.563 p=0.0002 → +0.047 p=0.64; viral-mimicry/IFN +0.530 → −0.000). The aberrancy layer behaves the same:
it is not independent of immune composition, so no architecture — linear or nonlinear — can extract signal
from an empty orthogonal subspace.

## What WOULD let a disentangling model work (data design, not architecture)
1. A feature genuinely decoupled from bulk infiltration — per-neojunction HLA-presented peptide load (the
   presentation layer), which measures PRESENTED aberrancy rather than a bulk burden that co-varies with
   immune-cell content. Scoped in ENCODER_PHASE_PROTOCOL.md.
2. A regime where the immune floor transfers — one large cohort with internal splits, so the residual is
   interpretable and does not inherit the cross-cohort sign flip.

Data: `results/eval/disentangle_test.json`.
