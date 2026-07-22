# Archive summary — evolutionary RNA-state / ICB

**Date:** 2026-07-22 · **Purpose:** single re-entry point if this project is ever reopened.
One page: what the hypothesis was, what was tested, the final verdict, what is dead, the one
path still open, and where every artifact lives.

---

## The hypothesis (verbatim, from NORTH_STAR)

> The captured bulk RNA-sequenced transcriptome of a human tumor contains sufficient information
> to reconstruct its latent evolutionary state in the context of immunotherapy response.

The mechanistic story: early driver mutations set a trajectory; genomic instability, epigenetic
remodeling, RNA-processing dysregulation and immune selection produce coordinated transcriptomic
abnormalities (splicing, intron retention, RNA editing, TE/ERV activation, fusions, cryptic
ORFs) that are downstream readouts of one latent evolutionary state **S**, which shapes
antigenicity and ICB response.

## Final verdict (2026-07-22)

**The strong form of the hypothesis is falsified on bulk RNA-seq.** Across many independent
attacks that do not share a failure mode, no tumor-intrinsic RNA representation carries ICB-
response information beyond immune composition, and the immune signal itself does not transport
across cohorts. The signal in bulk tumor RNA that tracks ICB response **is immune composition.**

This is a **rigorous cautionary negative with a mechanism attached**, not an inconclusive result.

## What is dead (tested, null, powered where noted)

| Layer | Test | Result |
|---|---|---|
| Antigen quantity (WES neoantigen proxy) | 416 samples / 5 cohorts | full-power null (perm p=0.78) |
| Non-reference RNA burden (splicing/IR/editing/TE/ERV) | within- and cross-cohort | no signal over immune floor; degrades it (multicollinear) |
| Frozen sequence encoders | EVA, HyenaDNA, Evo2-7B, **Caduceus** | **4/4 null** vs immune floor (Caduceus powered n=30, residual 0.343 perm p=0.86) |
| Trained shallow (MLP, fused blocks) | pooled + LOSO | null; nonlinear = linear |
| Trained deep VAE (scVI latent) | vs clonal ground truth, n≈60 | OOF ρ = permutation null, all 3 targets, both cohorts |
| dsRNA / viral-mimicry ΔG structure | IR-Alu, 286k pairs, expression-decoupled | **decoupled from IFN (novel) yet null** for response |
| Latent state vs clonality (literal claim) | purified malignant, WES-anchored, n≈60 | null (first properly powered test; prior was n=10) |
| Domain-invariant / robust learning (IRM, GroupDRO) | LOCO, the method class built for the sign-flip | ≤ ERM, ≤ chance — no invariant direction exists |

## The mechanism (why it fails — the real contribution)

The tumor-intrinsic "aberrancy" signal is a downstream readout of the inflamed-tumor state, and
its coupling to inflammation **sign-flips across cohorts** (GEP correlation +0.36 Gide →
−0.15 Hugo). That sign-flip is why the immune floor itself fails leave-one-cohort-out (canonical
3-cohort floor LOCO 0.59, perm p=0.044; ERM floor LOCO 0.399 — below chance, i.e. trained on two
cohorts it applies *backwards* on the third). Purity-corrected mediation: the phenotype→
infiltration path collapses to ~0 once purity is removed; only an antigen-presentation axis
retains a small tumor-intrinsic direct effect (log-odds +0.58, CI [+0.13,+1.15]) that does not
transport. Invariance methods confirm there is no cohort-stable tumor-intrinsic direction to
recover.

## The one path still open (substrate, not method)

Every null above purified the malignant compartment by **bulk deconvolution** (InstaPrism), an
estimate, not a measurement. The DL/SSL *method* space is now diligenced — frozen and trained,
supervised and unsupervised, capacity-scaling and invariance-constrained. What is left is not a
cleverer model but a cleaner substrate:

> **a matched single-cell / snRNA ICB melanoma cohort → build the tumor-intrinsic representation
> on directly measured malignant cells → then (optionally) the pre-registered fine-tune test
> (`EXPERIMENT_SCOPE_FINETUNE_RAWREADS_20260722.md`).**

Honest expectation: it confirms the negative (bulk composition dominates; single-cell removes
the deconvolution assumption but the biology likely holds). It is the only experiment that could
still flip the answer, and it requires new data + GPU, not more analysis on the current
substrate.

## Where everything lives

| Deliverable | Path |
|---|---|
| **This archive summary** | `docs/ARCHIVE_SUMMARY_20260722.md` |
| Gemini-exchange assessment + citation audit | `docs/GEMINI_ASSESSMENT_20260722.md` |
| dsRNA viral-mimicry + latent-state results | `docs/RESULTS_VIRALMIMICRY_LATENTSTATE_20260722.md` |
| Deep-learning diligence ledger (tested vs inference) | `docs/DEEP_LEARNING_DILIGENCE_20260722.md` |
| Pre-registered fine-tune experiment design | `docs/EXPERIMENT_SCOPE_FINETUNE_RAWREADS_20260722.md` |
| Rolling project status (full history) | `docs/PROJECT_STATUS.md` |
| Earlier manuscript + figure deck | `docs/WRITEUP.md`, `docs/manuscript.html`, `results/figure_deck.pdf` |
| Result JSONs | `results/eval/{dsrna_viral_mimicry_test,latent_state_clonality_probe,caduceus_frozen_test,invariance_loco_test}.json` |
| dsRNA feature map (checkpoint) | `results/features/iralu_pergene.parquet`, `dsrna_signature_persample.parquet` |
| Caduceus per-junction deltas | `results/features/caduceus_delta_gide{,32}.json` |
| Scorers | `analysis/score_iralu_dsrna.py`, `analysis/caduceus_score.py` |
| Two-prong figure | `results/fig_gemini_twoprong.png` |

## If reopening

1. Read this file, then `DEEP_LEARNING_DILIGENCE` (the method map) and `RESULTS_VIRALMIMICRY_LATENTSTATE` (the two sharpest tests).
2. The decision is binary: **write the negative-result + mechanism manuscript** (highest value, evidence is ready), **or** acquire single-cell data and take the one remaining shot.
3. Do not re-run bulk architecture search — it is diligenced. New substrate or a write-up; nothing in between adds information.
