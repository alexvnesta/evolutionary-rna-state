# Deep-learning diligence ledger

**Date:** 2026-07-22 · **Session:** `c55b6dde`
**Purpose:** an honest accounting of what deep learning has actually been *tested* on this
project vs. what remains *reasoned inference*. Written in response to the direct question
"did we do the full diligence there?" — the answer is **the frozen-embedding side is
thoroughly diligenced; the trained/fine-tuned side is a judgment call, not a demonstrated
null.** This ledger makes that boundary explicit so it is documented, not implicit.

The classification below is deliberately conservative: a row is **TESTED** only if a model was
run on project data against the immune-floor positive control (LOCO 0.767) or a permutation
null; it is **INFERENCE** if the verdict rests on an argument, however well-grounded, rather
than a run.

---

## A. TESTED — run on project data, null (high confidence)

| # | Model | Class | Input | Result | Source |
|---|---|---|---|---|---|
| 1 | **EVA** (1.4B RNA LM) | frozen embedding | bulk expression | Collapsed to linear function of expression, R²=1.000000; LOCO 0.413 vs PCA 0.404 | `EVO2_ENCODER_RESULT`, `FINAL_VERDICT` |
| 2 | **HyenaDNA** | frozen embedding | sequence | LOCO 0.421 (chance; perm p=0.77) | `FINAL_VERDICT` |
| 3 | **Evo2-7B** | frozen, zero-shot likelihood | novel-junction seqs (1024 bp, spliced vs contig) | Within-Gide n=32 residual 0.354 (perm p=0.83); within-Hugo n=22 residual 0.576 (p=0.259) | `EVO2_INTERACTION_INSIGHT`, `evo2_two_block_gide32.json` |
| 4 | **Supervised MLP** | trained (shallow) | fused abundance(64) + aberration(35) | Null and *worse than linear*: pooled 0.61→0.46 adding aberration; MLP no better than linear | Gemini-session run, `GEMINI_ASSESSMENT` |
| 5 | **scVI** (deep VAE) | trained, unsupervised | expression | Latent = permutation null vs clonal ground truth, both cohorts (perm p 0.49–0.93) | `RESULTS_VIRALMIMICRY_LATENTSTATE`, this session |

**Two recurring failure modes:** (i) **collapse to expression** (EVA, R²=1.0); (ii) **track
composition / carry nothing tumor-intrinsic** (HyenaDNA, Evo2, scVI). Against the immune-floor
positive control (LOCO 0.767), every one of these is at chance. Note that #5 is a *deep,
trained, unsupervised* model — so "we only tested frozen models" is not accurate; scVI was
trained end-to-end and still null against a tumor-intrinsic anchor.

## B. TESTED — run this session, closes a named gap

| # | Model | Class | Input | Result | Source |
|---|---|---|---|---|---|
| 6 | **Caduceus-1k** (`kuleshov-group/caduceus-ps_seqlen-1k`) | frozen embedding-delta | Gide-n32 junction seqs (same substrate as Evo2 #3) | Ran on Modal (H100 job `9790fc35` + A100 job `c7e1f54b`). **Powered n=30:** floor 0.852, Caduceus alone 0.684, floor+Caduceus 0.826 (Δ **−0.026**); residual-on-floor AUROC **0.343, perm p=0.860**. Essentially identical to Evo2 within-Gide (residual 0.354, perm p=0.83). | `caduceus_frozen_test.json` |

> **#6 is powered and complete.** The n=30 test uses the same `gide32_top200` / 1864-junction substrate as the Evo2 within-Gide test, so it is a direct apples-to-apples comparison — and it lands in the same place: frozen Caduceus *subtracts* from the immune floor and its residual is a dead-center null (perm p=0.86). This closes the frozen-encoder sweep at **4-for-4 null, all powered** (EVA, HyenaDNA, Evo2, Caduceus) against the immune floor. (Modal H100/A100 reserve was flaky during this session — several submits timed out at the reserve stage before two got through; the two that ran gave identical direction, n=11 and the powered n=30.)

## C. INFERENCE — not run; verdict rests on argument (lower confidence, honestly flagged)

| # | Approach | Why deferred | The argument for expecting null | Confidence |
|---|---|---|---|---|
| 7 | **PEFT / fine-tuning** a 1–7B sequence LM end-to-end | Gated in `ENCODER_PHASE_PROTOCOL` on n≈100 being too small to fine-tune without overfitting | Any model trained *toward response* optimizes toward the composition confound that lives in the label; the MLP test (#4) confirmed more capacity doesn't help — nonlinear = linear = null | medium |
| 8 | **Masked / pathway-informed VAE** (P-NET / bMAE direction) | Assessed in `GEMINI_ASSESSMENT`, not built | On composition-dominated bulk, a pathway-masked latent re-learns the confound (the "viral mimicry" node covaries with the immune node because that is the true biology); same substrate as the null scVI test | medium |
| 9 | **End-to-end training on raw reads** (the original ambition) | Never run at cohort scale; only frozen embeddings of pretrained raw-read models tested | Same label-confound logic as #7; and frozen Evo2 already collapsed | medium-low |

## D. TESTED — the invariance/robustness class, aimed directly at the sign-flip (run this session)

The project's wall is **non-transportability**: the tumor-intrinsic signal's coupling to
response sign-flips across cohorts, and the immune floor itself fails LOCO (canonical 3-cohort
floor LOCO 0.59, perm p=0.044 — "the floor itself fails LOCO", `PROJECT_STATUS`). Buckets A–C
either ignore that (fit harder / bigger) or re-learn the confound. There is one method class
built *specifically* for train-on-some-cohorts-transport-to-another: **domain-invariant /
distributionally-robust representation learning.** It was never tried until now.

| # | Method | Principle | Result (LOCO mean, n=65, 3 cohorts) | 
|---|---|---|---|
| 10 | **ERM** (baseline) | pooled logistic | floor 0.399 · scVI latent 0.442 |
| 11 | **IRM** (IRMv1) | penalize env-varying optimal classifier | floor 0.328 · scVI 0.422 |
| 12 | **GroupDRO** | minimize worst-environment loss | floor 0.330 · scVI 0.475 |

**Result: the invariance class does not rescue transportability.** IRM and GroupDRO are ≤ ERM
and ≤ chance on both the immune floor and the scVI RNA latent. This is the expected outcome
under the mechanism — if the signal sign-flips because it is a downstream readout with
cohort-specific coupling, then *there is no invariant direction to find*, and a correct
invariance method drops it rather than inventing one. The ERM floor LOCO of 0.399 (below chance)
is the sign-flip signature shown directly: a model trained on two cohorts applies backwards on
the third. `invariance_loco_test.json`. **This closes Bucket 3** — the one refinement class that
targeted our actual failure mode, rather than variance or capacity, also finds nothing.

---

## The genuinely untested frontier

Collapsing A–C: the one combination no test has covered is

> **trained (not frozen) × raw-read model × tumor-intrinsic objective (not response) × on a
> purified or single-cell malignant compartment (not bulk).**

Every existing test violates at least one of those conditions:
- #1–3, #6 are *frozen*, not trained.
- #4, #7–9 train *toward response* (label-confounded).
- #5 trains unsupervised but on *bulk* expression (composition-dominated) and reconstructs
  expression, not a tumor-intrinsic evolutionary target.
- #10–12 (invariance class) operate on bulk representations and correctly find no invariant
  direction — they close the *method* gap (we tried the class built for our failure mode) but
  not the *substrate* gap (still bulk, not purified/single-cell).

Note the frontier is now defined by **substrate, not method**: we have tried frozen and trained,
supervised and unsupervised, capacity-scaling and invariance-constrained. What is left is not a
cleverer model — it is running any of these on a *directly measured* malignant compartment
instead of bulk. That is why the remaining path requires single-cell data, not more architecture
search.

**My expectation** is that this frontier also collapses — frozen Evo2 already collapsed, and
bulk composition dominates every representation we have built. But this is **inference, not a
demonstrated null**, and I will not present it as equivalent to rows A. It is the honest edge
of the project's evidence.

## What would actually close it

1. **Cheap partial closure (this session):** Caduceus *frozen* embeddings on the Evo2 junction
   set — row #6. Completes the frozen-encoder sweep. Likely null, but a real data point.
2. **The real test (needs new substrate + GPU):** fine-tune a raw-read model (Evo2 / Caduceus)
   with a *reconstruction / tumor-intrinsic-phenotype* objective on **single-cell or
   deconvolution-purified malignant** input — never toward response. This is the natural
   companion to the single-cell path (see `RESULTS_VIRALMIMICRY_LATENTSTATE` §"what this does
   not settle"), because both require the malignant compartment resolved properly rather than
   estimated from bulk. It is the experiment most likely to *either* change the answer *or*
   overfit spectacularly at n≈100 — which is exactly why it needs single-cell scale to be
   interpretable.

**Bottom line.** Deep learning is not un-diligenced — five model classes including a trained
deep VAE are demonstrated null against the right controls. But the specific cell
[trained × raw-read × tumor-intrinsic × purified] is open, and the verdict there is currently
an argument, not a result. Closing it fully requires single-cell data and GPU, not more
modeling on the current bulk substrate.
