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
- The residualization was fold-contained throughout (floor→evo regression fit on TRAIN only inside
  each CV split). The eye-catching 0.643 at n=13 was a **single-seed (seed0)** value of that same
  procedure; its **20-seed mean was already 0.456** (chance), with permutation p=0.190. In other words,
  0.643 was single-split variance at tiny n, not a data-leakage artifact — averaging over seeds and
  doubling n to 32 (20-seed mean 0.354, p=0.83) both confirm the null. Report the multi-seed mean, not a
  single lucky split, at small n.
- Two H100 scoring jobs were orphaned by submit-helper timeouts (Friday queue congestion,
  not a code fault; actual scoring is ~2-3 min). Switching to A100-80GB resolved it.

## Independent replication — within-Hugo (n=22, 12R/10N)
Ran the identical pipeline within a second cohort (Hugo 2016), deepened to 22 samples (from 17) as more
Hugo junction extractions became available. Scored the Hugo-specific novel junctions on A100.

| Block | AUROC (20-seed mean) |
|---|---|
| Immune floor | 0.593 |
| Evo2 alone | 0.502 (chance) |
| Floor + Evo2 | 0.593 (no gain, Δ 0.000) |
| Evo2 residualized on floor (fold-contained) | 0.576 |
| — permutation p (residual) | **0.259 (not significant)** |

**Same null verdict in a second cohort — confirmed by permutation.** The residual reads 0.576 (20-seed mean;
0.586 at the 10-seed setting used for the permutation), but this is NOT significant: the cohort-internal
permutation null gives p=0.259 (observed 0.586 vs null mean 0.494), and the per-seed
residual ranges 0.458–0.700 (std 0.065) — i.e. a single CV seed can land anywhere from chance to 0.70 at
n=22. Critically, the direct tests agree it is null: Evo2 alone is exactly chance (0.502) and adding Evo2 to
the floor yields zero gain (0.593→0.593). This is the same small-sample-variance signature the anti-collapse
control was built to catch (cf. the n=13 Gide 0.643 that vanished on seed-averaging + more n). Hugo's own
floor is weak (0.593) — consistent with the Phase 3 finding that Hugo's response biology differs — but the
permutation test controls for that, and Evo2 still carries no significant independent signal. Two independent
within-cohort tests (Gide n=32 p=0.83, Hugo n=22 p=0.259) now agree: the Evo2
novel-junction aberrancy layer carries no independent ICB-response signal.
