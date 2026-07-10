# Final verdict: "Evolutionary RNA state → antigenicity → ICB response"

Consolidated across the full investigation. This document is the top-level read; the
per-thread notes (io_multicohort_clonality.md, io_move2_autorun_watch.md,
SYNTHESIS_hypothesis_reconsideration.md) hold the detail. Figure: fig_master_verdict.png.

> **⚠️ CORRECTION NOTICE (2026-07-09) — a four-track forensic re-audit found several verdicts below
> are not supported by what was actually tested. See `AUDIT_SYNTHESIS.md`. Key corrections:**
> 1. **The "LEARNED REPRESENTATION" row is not a foundation model.** It is PCA on the top-2000 variable
>    genes — a linear summary of bulk expression. The "LOCO AUROC 0.49" attributed to Orthrus/encoders is
>    this PCA/scVI baseline; **Orthrus was never run at patient scale** (only a 5-gene timing benchmark exists).
>    The EVA "encoder" feature is provably a linear function of expression (R²=1.0 from 39 expression PCs).
>    No foundation model was ever fed a patient-specific or aberrant sequence. → **verdict should read "never validly tested"**.
> 2. **The 416/5-cohort "RNA antigen QUANTITY" row is a WES (DNA) exome-neoantigen proxy on ≤264 samples**,
>    not direct transcriptome splice/ERV/fusion counts. The DEAD label is valid *for the DNA proxy*; it does
>    not test the RNA phenotypes. The non-reference RNA features were quantified on **0/106** samples.
> 3. **The DEAD verdicts at n=40 / n=66 are underpowered** (≈39% power for true AUROC 0.65; learned_rep
>    95% CI [0.385, 0.743] spans chance→useful) — "no signal detectable at low power," not "no signal exists."
> 4. **The immune floor was regressed out as a confounder without testing whether it is a confounder or a
>    mediator.** A first mediation analysis (n=40, marker-proxy mediator) suggested a mediated chain; but a
>    corrected re-analysis (2026-07-10, n=106, using a REAL reference-based deconvolution mediator — InstaPrism
>    on the Pozniak 2024 melanoma atlas with a malignant compartment) supersedes it. With a real mediator and
>    tumor-purity correction: (a) the phenotype→infiltration a-path is largely a **purity artifact** — it
>    collapses to ~0 for both axes once purity is removed; (b) the n=40 viral-mimicry/IFN "mediation" does NOT
>    survive; (c) **but the antigen-presentation axis retains a tumor-intrinsic DIRECT effect on response**
>    (log-odds +0.58, 95% CI [+0.13, +1.15], P>0=0.995), independent of purity and infiltration. So "adjust for
>    immune, nothing survives" is wrong as originally stated: one tumor-intrinsic expression axis does survive
>    a proper purity+infiltration adjustment. It is still an expression-derived axis (cannot beat the immune
>    floor as a predictor; value is mechanistic). See `AUDIT_SYNTHESIS.md` §Phase 1b and
>    `mediation_purity_corrected_n106.png`.

## The hypothesis (as stated by the author)
"The captured bulk RNA-sequenced transcriptome of a human tumor contains sufficient
information to reconstruct its latent evolutionary state in the context of immunotherapy
response" — via a HYBRID of learned representations from raw RNA-seq reads + interpretable
transcriptomic features (expression, splicing, RNA editing, TE activity, fusions, cryptic ORFs).

Embedded commitments: (C1) the RNA abnormalities are COORDINATED (one latent state);
(C2) that state drives ANTIGENICITY and thereby ICB RESPONSE; (C3) it can be RECONSTRUCTED
from bulk RNA, including via a LEARNED REPRESENTATION.

## Scorecard
| Operationalization | Data | Verdict | Key statistic |
|---|---|---|---|
| RNA antigen QUANTITY (WES exome-neoantigen proxy, *not* RNA counts) | ≤264, WES | **DEAD (for DNA proxy)** | PC1 at null p=0.77; LOCO AUROC 0.44-0.50 < TMB 0.62 — see correction note |
| RNA load × CLONALITY interaction | 189, 3 coh | **DEAD** | pooled LR p=0.82; n=40 Riaz hint was noise |
| Tumor-intrinsic RNA over IMMUNE floor (IFN) | 40, 2 coh | **DEAD** | joint p=0.37; 2 feats lingered p~0.10 |
| Tumor-intrinsic over HARDENED floor (11 cell types) | 40, 2 coh | **DEAD** | joint p=0.54; lingering feats absorbed to p=0.28-0.34 |
| "LEARNED REPRESENTATION" (= expression PCA, NOT a foundation model) vs immune floor | 66, 3 coh | **NEVER VALIDLY TESTED** | LOCO AUROC 0.49 is the expr-PCA/scVI baseline (95% CI [0.35,0.63], underpowered); "degrades floor" is a multicollinearity artifact (92% collinear); Orthrus never run — see correction note |
| — immune floor (positive control) | 66, 3 coh | **WEAK** | LOCO AUROC 0.59, perm p=0.044 (was 0.70 at 2 coh) |
| Heterogeneity (MATH) → response | 189, 3 coh | **REAL** | p=0.014 — but DNA-side, not RNA |
\*n=40 (2 of 5 cohorts); the powered re-test fires automatically when liu2019 lands (~n=132).

## What is definitively established (full-power, not provisional)
1. RNA antigen QUANTITY carries no ICB-predictive signal. Three independent operationalizations
   (marginal counts, count×clonality, count-over-TMB) on the full 416-sample / 5-cohort data
   are null. C2 fails for antigen counts. This is the strongest, most durable result.
2. The ICB signal in bulk RNA IS immune composition. A multi-cell-type deconvolution block
   (11 xCell-style populations) predicts response through 8/11 immune populations (checkpoint
   ρ=+0.52, CD8 T ρ=+0.47, NK/Treg/cytotoxic significant), while non-immune compartments
   (stroma, endothelial) are flat. After removing this floor, NO tumor-intrinsic RNA feature
   retains signal (joint p=0.54 at n=40; p=0.383 at n=66/3 cohorts). At 3 cohorts with proper
   cohort-adjustment the surviving immune populations are the T-cell/checkpoint axis specifically
   (checkpoint p=0.009, CD8 T p=0.028, Treg p=0.031, CD4 p=0.038; 4/11 significant, stroma/endothelial
   still non-significant and negative) — the signal is immune-specific, not generic infiltration.
   This confirms, at this scale, the CIBERSORTx/xCell critique
   raised by a public commenter (E. Garcia Lecaros) on the author's own LinkedIn post announcing
   the project — the exact objection that bulk-RNA ICB signal reflects inferrable immune-cell
   composition rather than reconstructed tumor evolutionary state.

## What was newly tested this run — the hypothesis's OWN stated method (C3)
The learned-representation half had never been tested; every prior analysis used hand-engineered
features only. Built unsupervised RNA embeddings (PCA/kernel-PCA/MLP-autoencoder, fit inside CV
to prevent leakage) and benchmarked against the immune floor under two schemes:
  - Random 5-fold CV: learned rep 0.75 ≈ immune floor 0.76 — looks competitive.
  - Leave-one-cohort-out (the honest test): learned rep collapses to 0.56 (perm p=0.153, does NOT
    beat chance across cohorts) while the immune floor holds at 0.70 (perm p=0.007). Adding the
    representation to the floor DEGRADES it (0.70 → 0.56). [All perm p from headtohead_repr_40samp.json;
    500-perm null. Values are stable to ±0.01 across Monte Carlo draws — both sides of p=0.05 are unambiguous.]
THREE-COHORT UPDATE (hugo2016 expression landed mid-run; auto-runner fired, hugo crosswalk built
and validated 26/28 → iAtlas): with a genuine third held-out cohort the picture sharpens. Learned
rep LOCO AUROC drops to 0.49 — exactly chance (perm p=0.41); the immune floor's cross-cohort signal,
tested honestly against a truly independent third cohort, shrinks to 0.59 (perm p=0.044) from the
optimistic 2-cohort estimate of 0.70. floor+rep = 0.50 (representation drags the floor to chance).
So the 2-cohort head-to-head was if anything GENEROUS to the hypothesis; adding hugo confirms the
learned representation carries no transferable signal and the real immune-composition effect is
modest under strict validation. (headtohead_repr_66samp.json.)

The random-CV parity was within-cohort overfitting (learned-rep score variance 0.060 vs floor
0.010, 6×). ROBUSTNESS: a 9-config hyperparameter sweep (5/10/20 PCs × 1000/2000/5000 genes;
learned_rep_sensitivity_40samp.json) holds the learned-rep LOCO AUROC in 0.56-0.63 for EVERY
config — never approaching the immune floor's 0.70 — while random-CV stays 0.69-0.79 throughout,
confirming the overfitting gap is systematic, not a lucky configuration. VERDICT on C3 at current
power: a learned RNA representation carries no cross-cohort ICB signal beyond immune composition,
and actively harms the predictor that works. This is the first honest test of the method the
hypothesis actually proposed, and it is robust to representation hyperparameters.

## What is NOT yet killed (the honest remaining gaps)
1. POWER. hugo2016 expression landed mid-run, taking the head-to-head to n=66 / 3 cohorts (gide,
   hugo, riaz) — already a THIRD held-out cohort for LOCO, and it strengthened the null. liu2019
   (n=122, has TMB) remains the big step: it takes the response tests to n~180 and a FOURTH held-out
   cohort. The extended auto-runner (move2_autorun.py v6, now with a validated hugo crosswalk) fires
   ALL tests automatically when liu lands; a monitor is polling for it.
2. THE LITERAL LATENT VARIABLE. The hypothesis says "reconstruct the evolutionary STATE"
   (clonal architecture), not "predict response." Direct probe RNA-axis vs clonality is capped
   at n=10 (RNA∩clonality, riaz-only) until liu expression lands; at n=10 no RNA axis tracks
   clonality (all |ρ|<0.19) but that is underpowered. liu unblocks this to n~132 too.
3. REPRESENTATION CAPACITY. The learned rep here is deliberately low-capacity (PCA/shallow AE)
   because n=40 forbids anything deeper without memorizing. A raw-read foundation-model embedding
   (the author's literal "learned representations from raw reads") could differ — but only a
   larger sample makes that testable without overfitting. At n=40 the low-capacity result is the
   scientifically defensible one.

## Bottom line
The hypothesis is NOT vindicated and its antigen-quantity core is firmly dead. Its stated
learned-representation method, tested honestly for the first time, fails at cross-cohort transfer
and adds nothing over immune deconvolution. What survives is narrow and orthogonal to the original
claim: a coordinated, proliferation-independent RNA-regulator activity state that MARKS immune-cold
tumors (does not drive antigenicity), and a DNA-side heterogeneity effect on response. The single
scientifically live question — whether ANY of this reaches significance at proper power — is fully
wired to auto-answer when liu2019 expression lands. If the powered tests come back null (the
current n=40 trajectory strongly predicts they will), the hypothesis is dead in full.

RECOMMENDATION: write up the antigen-quantity negative and the immune-composition-dominance result
now (both full-power). Hold the RNA-state re-scope as provisional pending liu. Do not invest in
higher-capacity learned representations until n is large enough to test them without overfitting.

## Artifacts (this autonomous run)
- fig_master_verdict.png — 3-panel: head-to-head, immune deconvolution, scorecard
- headtohead_repr.py / headtohead_repr_40samp.json — the learned-rep vs floor benchmark (LOCO perm)
- learned_representations_40samp.csv — the PCA/kPCA/AE embeddings
- immune_deconvolution_block_40samp.csv — 11-cell-type xCell-style scores
- residual_over_hardened_floor_40samp.json — tumor-intrinsic over multi-cell-type floor
- move2_autorun.py v5 — extended runner (activity + reframing P1-P3 + head-to-head + ID diagnostics),
  wired for liu via passthrough crosswalk; move2_watch_live.py monitoring
Upstream: SYNTHESIS_hypothesis_reconsideration.md, io_multicohort_clonality.md, io_move2_autorun_watch.md
