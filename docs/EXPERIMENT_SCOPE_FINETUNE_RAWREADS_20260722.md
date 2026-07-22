# Experiment scope — fine-tuning a raw-read model toward the tumor-intrinsic state

**Date:** 2026-07-22 · **Session:** `c55b6dde` · **Status:** design / scoping (not yet run)
**Frontier it targets:** the one untested cell in `DEEP_LEARNING_DILIGENCE_20260722.md` —
*trained (not frozen) × raw-read model × tumor-intrinsic objective (not response) × purified/
single-cell compartment.* This document specifies the experiment concretely enough to build,
including — critically — **what result would falsify the hypothesis and what would support it.**

---

## 1. The design in one paragraph

Take a raw-read genomic LM (Evo2-7B, weights already hydrated on Modal; or Caduceus if the
frozen test motivates it), and instead of freezing it or training it toward ICB response,
**fine-tune it (PEFT/LoRA) to reconstruct a tumor-intrinsic evolutionary target** — clonal
architecture (heterogeneity_index / subclonal_fraction) or the malignant-cell expression
program — from the sample's aberrant-transcript sequences. Response is never in the loss. The
test is then the same two-block protocol as everything else: does the fine-tuned representation
carry information about ICB response *over the immune floor*, and does it transport across
cohorts (LOCO)?

## 2. Why this is the right next DL experiment (and the only one worth GPU)

- It is the **only** design that avoids the label-confound that killed the MLP (#4) and gates
  the PEFT/VAE rows (#7–8): the model never sees response, so it cannot rediscover composition
  through the label.
- It is the **only** design that trains rather than freezes on raw reads — closing the gap that
  frozen Evo2/HyenaDNA/Caduceus leave open.
- It couples to the single-cell path: the tumor-intrinsic target is only clean on a resolved
  malignant compartment, so the two efforts share a substrate.

## 3. Inputs required (and their status)

| Input | Status | Gap |
|---|---|---|
| Aberrant-transcript sequences per sample (junctions, retained introns, TE loci) | **On disk** for junctions (`results/junctions/junction_seqs*.json`, 1024 bp windows) | IR / TE-locus seqs need extraction from the 106 STAR BAMs |
| Tumor-intrinsic target: clonal architecture | **On disk**, within-cohort (Hugo 27 CCF, Riaz 33 VAF-proxy) | Not cross-cohort comparable; Gide has no WES |
| Tumor-intrinsic target: malignant expression program | InstaPrism Malignant fraction on disk (bulk estimate) | **Single-cell would be far cleaner** — bulk deconvolution is the weak link |
| Raw-read model weights | Evo2-7B **hydrated on Modal** (`claude-science-evo2-weights`, H100) | Caduceus weights not yet on Modal; mamba-ssm CUDA image not built |
| GPU | Modal H100 available (`byoc:modal`); box GPU held by PubMedBERT service | ~$22 of $30 Modal budget left — enough for a scoped LoRA run, not a large sweep |

## 4. The protocol

1. **Target definition.** Primary target = within-cohort clonal architecture (the literal
   evolutionary state). Secondary = malignant-cell expression program (needs single-cell to be
   credible). Both are tumor-intrinsic; neither is response.
2. **Fine-tune.** LoRA adapters on Evo2-7B, objective = regress the sample's pooled
   aberrant-sequence representation onto the target. **Fold-contained**: adapters trained on
   train-fold samples only, target standardized in-fold. Cohort-internal to respect the
   non-comparability of Hugo/Riaz clonality.
3. **Extract** the fine-tuned per-sample representation (OOF on held-out fold).
4. **Two-block test** (the project-standard protocol): floor / representation-alone /
   floor+representation, 20-seed CV within cohort, cohort-internal permutation null, and LOCO
   for transportability.

## 5. Pre-registered falsification criteria (decide BEFORE running)

This is the part that makes it worth doing — a pre-committed decision rule so the result is
interpretable either way and not a garden of forking paths:

- **SUPPORTS the hypothesis** if, and only if: the fine-tuned representation (a) reconstructs
  the tumor-intrinsic target above its permutation null within cohort (OOF ρ clearly off the
  null band, unlike scVI's ρ=null in Prong B), **AND** (b) adds AUROC over the immune floor
  with a lower CI bound above 0, **AND** (c) that increment does not vanish under LOCO. All
  three. Any one failing = not support.
- **FALSIFIES (strengthens the negative)** if: the representation is null against its own
  tumor-intrinsic target (like scVI), OR it reconstructs the target but adds nothing over the
  floor, OR it adds within-cohort but collapses under LOCO (the sign-flip signature).
- **Expected outcome (stated honestly, in advance):** falsification — most likely the
  representation reconstructs the target weakly and still adds nothing transportable, because
  bulk composition dominates and n≈100 invites overfitting. The value is that this is the
  *last* DL configuration, so a null here converts "we expect DL fails" into "we showed DL
  fails, trained and frozen, supervised and unsupervised, on the best available objective."

## 6. Cost, risk, and the honest recommendation

- **Cost:** one scoped LoRA run on Modal H100 ≈ a few GPU-hours, within remaining budget. Image
  build (mamba-ssm/PEFT) is the main setup cost.
- **Risk:** at n≈100 bulk, a within-cohort "support" result is more likely overfitting than
  signal — which is exactly why criterion (c), LOCO transportability, is non-negotiable.
- **Recommendation:** **do not run this on bulk alone.** Its interpretable version needs the
  single-cell malignant compartment. The right sequence is: (1) close the cheap frozen-Caduceus
  gap now; (2) acquire a matched single-cell ICB cohort; (3) *then* run this fine-tune with the
  malignant program as target. Running it on bulk now would produce a likely-null or
  likely-overfit result that neither supports nor cleanly falsifies — the worst kind of
  experiment. Documented here so the design is ready when the single-cell substrate exists.
