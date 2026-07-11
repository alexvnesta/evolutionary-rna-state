# Evo2 Novel-Junction Encoder Pass — Within-Gide Result

**Date:** 2026-07-11  **Status:** COMPLETE (negative)  **Compute:** Modal A100/H100, Evo2-7B bf16

## Question
Does a genomic foundation model (Evo2-7B) scoring the *sequence aberrancy* of novel
splice junctions add ICB-response signal beyond the reference immune floor?

## Feature
Per sample, the top-200 novel splice junctions by read support (canonical chroms,
≥10 reads, absent from GENCODE v46 introns [convention A, 72.8% known-match validated],
recurrent in ≥2 Gide samples). Each junction scored by Evo2-7B as
delta = mean-LL(spliced donor+acceptor window) − mean-LL(contiguous reference window).
Per-sample block = {mean, median, min, frac_neg, read-weighted-mean} of its junction deltas.
1,887 unique junctions scored across 32 samples.

## Evaluation frame
Within-Gide only (n=32, 17R/15N). Cross-cohort LOCO is invalid for this 3-cohort set:
the immune floor itself fails to transfer (Phase 2, LOCO 0.507), so no block is
interpretable cross-cohort. Within-Gide is the only frame with a working positive control.

## Result (20-seed mean 5-fold CV AUROC)
| Block | AUROC |
|---|---|
| Immune floor (positive control) | **0.792** |
| Evo2 junction-aberrancy alone | 0.513 (chance) |
| Floor + Evo2 | 0.768 (Δ −0.024, degrades) |
| Evo2 residualized on floor (fold-contained) | 0.354 |
| Permutation p (residual) | 0.83 |

## Verdict
The Evo2 novel-junction aberrancy layer carries **no independent ICB-response signal**
within Gide at n=32. Adding it degrades the immune floor. This is consistent with every
prior non-reference layer tested (editing, IR, TE, splicing, learned expression encoders):
the immune-composition floor is not beaten by any non-reference RNA feature at this scale.

## Rigor notes
- An in-sample (non-fold-contained) residualization gave a spurious 0.643 at n=13;
  fold-contained residualization and doubling n to 32 both collapse it to chance —
  a textbook small-sample optimism artifact, caught by the pre-registered anti-collapse control.
- Two H100 scoring jobs were orphaned by submit-helper timeouts (Friday queue congestion,
  not a code fault; actual scoring is ~2-3 min). Switching to A100-80GB resolved it.
