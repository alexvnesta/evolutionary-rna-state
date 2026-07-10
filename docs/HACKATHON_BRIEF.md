# Hackathon result brief — the first honest test of the non-reference RNA layer

**Project:** evolutionary-rna-state · **Session:** `23cf8106` (OPERON) · **Date:** 2026-07-10
**Hackathon:** Built with Claude: Life Sciences (Jul 7–13, 2026)

## The one-sentence claim
We built the project's first **cohort-scale non-reference RNA feature matrix** and ran a **pre-registered
two-block test** of whether those features — RNA editing, intron retention, splicing burden, TE/ERV family
activity — carry immune-checkpoint-blockade (ICB) response signal **beyond** an established immune-composition
floor. This is the first time the hypothesis's *own* features have been tested; every prior verdict in this
project used expression, expression-PCA, or a DNA antigen-quantity proxy — never the non-reference layer.

## Why this is not just COMPASS
COMPASS (Shen 2025, *Nat Med*, DOI 10.1038/s41591-026-04502-7) is a concept-bottleneck transformer over
**reference gene expression** (44 immune concepts, 10,184 tumors). Its signal is immune-composition-centric.
Our hypothesis is orthogonal: that **non-reference** transcriptome events — the readout of genomic instability,
RNA-processing dysregulation, and TE de-repression — encode a latent *evolutionary* state that expression
alone misses. Testing that requires features COMPASS never computes. Whether or not the signal is there, this
is the question COMPASS leaves open, not a re-run of it.

## What was built (all committed this session)
1. **Unified Apple-Silicon Nextflow workflow** (`pipelines/main.nf`, `-profile apple_silicon`): one DSL2
   entry that fans a spine BAM channel to the validated editing / intron-retention / TE-ERV subworkflows.
   Compiles clean on Nextflow 26.04. (`README_UNIFIED.md`)
2. **First cohort-scale non-reference matrix** (`results/predictor/nonref_matrix_cohort.parquet`): 25
   samples × 66 features (editing 4, IR 4, splice 1, TE/ERV family 57), across 3 cohorts (Gide 14, Hugo 6, Riaz 5) — deepened as Hugo alignment landed in the background.
   Built by running the four callers on existing HISAT2 BAMs (`callers_on_bam.sh` / `fast_callers.sh`).
3. **Pre-registered evaluation** (`docs/EVAL_PROTOCOL.md` + `analysis/two_block_eval.py`), locked with a
   fixed success criterion, CV frame, and four rigor checks BEFORE any AUROC was computed.

## Result
See `results/eval/nonref_vs_floor_grouped5fold.json`, `rigor_checks.json`, `two_block_result.png`.

| block | features | AUROC | 95% CI | perm p |
|---|---|---|---|---|
| A — immune floor | 5 | 0.792 | [0.54, 1.00] | 0.021 |
| B — non-reference | 66 | 0.333 | [0.04, 0.62] | 0.663 |
| C — floor + non-ref | 71 | 0.542 | [0.23, 0.86] | 0.299 |

**ΔAUROC (C − A) = -0.250.** CV frame: within-Gide 5-fold (grouped by patient); Riaz n=5 & Hugo n=6 below the pre-set ≥10 hold-out threshold. n = 14 (6R/8N). (Full matrix n=25; primary frame is Gide-only, so its n is unchanged — the 4 new Hugo samples feed the 3-cohort LOCO below.)

## Honest verdict — NEGATIVE, and clean
**The non-reference RNA layer carries no independent ICB-response signal at this n, and adding it degrades
the immune floor.** This is consistent across all three analyses:

1. **Primary (within-Gide, n=14):** floor AUROC **0.792** (perm p=0.021 — the established predictor
   replicates); non-ref alone **0.333** (p=0.663, indistinguishable from chance); combined **0.542**,
   **ΔAUROC(C−A) = −0.250** — non-ref *dilutes* the floor.
2. **Composition check (purity covariate):** adding InstaPrism tumor-purity barely moves anything
   (floor+nonref 0.479 → +purity 0.458; non-ref alone 0.313 with or without purity). The non-ref failure
   is **not** a purity confound — it is a genuine absence of signal.
3. **Secondary 3-cohort LOCO (Gide/Hugo/Riaz, n=25):** floor LOCO **0.571**, non-ref LOCO **0.474**,
   floor+non-ref **0.583** — non-ref at chance and adding nothing meaningful even in a genuine
   leave-one-cohort-out frame across all three cohorts. (The floor's LOCO drops from its within-cohort 0.79
   because cross-cohort transfer at n=25 is hard — but non-ref never contributes.)

**Against the pre-registered success criterion** (B beats its null AND ΔC−A CI excludes 0): **not met on
either count** — a decisive negative.

**Why this is a real result, not a null-by-underpowering artifact:** the *positive control* (immune floor)
is significant at the very same n (p=0.021 within-Gide). A frame that recovers the known
predictor but not the non-ref block is informative about the non-ref block, not merely about power.

**Context:** this matches the honest prior — the ADAR editing regulon LOCO collapsed 0.647→0.535 (p=0.31)
at n=106, and the immune floor was already the real cross-cohort predictor (LOCO 0.767, p=0.001). But
unlike every prior "DEAD" verdict in this project (which used expression / expression-PCA / a DNA
antigen-quantity proxy), **this is the first test run on the hypothesis's actual non-reference features.**
It is a legitimate, reportable negative — bounded by n — not a repeat of the proxy-based verdicts.

## Scope and limits (stated plainly)
- **n is small** (14 samples). This is a first matrix built under a 3-day deadline on a single Apple Silicon
  machine, not the full n=106. The result bounds the effect size detectable at this n; it does not settle the
  hypothesis.
- **The learned-representation (encoder) half is NOT tested here** — scoped as future work
  (`docs/ENCODER_PHASE_PROTOCOL.md`) precisely because a rushed run would reproduce the `Wn @ E`
  expression-collapse (R²=1.000) the forensic audit caught.
- Locus-level TE (Telescope) is family-level here; escalation to loci is a documented cloud step.

## Reproduce
```bash
source pipelines/env.sh
# 1. features (existing BAMs):
bash results/nonref_run/fast_callers.sh <acc> <cohort> <bam>   # per sample; or the batch loop
python analysis/assemble_nonref_matrix.py
# 2. pre-registered test:
python analysis/two_block_eval.py --nonref results/predictor/nonref_matrix_cohort.parquet --out results/eval
python analysis/plot_two_block.py results/eval/nonref_vs_floor_*.json results/eval/two_block_result.png
```
