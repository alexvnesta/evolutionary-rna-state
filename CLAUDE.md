# CLAUDE.md — cross-agent working notes

Read this first. It tells any agent (or future session) what is live, what is
settled, and where the load-bearing caveats are. Keep it short and current;
put detail in `docs/`.

## The one caveat that governs everything

The core hypothesis has already been tested on its own terms and **does not
survive at cross-cohort transfer** (see `docs/FINAL_VERDICT.md`):

- A learned tumor-intrinsic RNA representation is at **chance** across cohorts
  (leave-one-cohort-out AUROC **0.49**, perm p=0.41) and **degrades** an
  immune-composition floor when added to it.
- The transferable ICB signal in bulk melanoma RNA **is immune composition**
  (T-cell / checkpoint axis), not reconstructed tumor evolutionary state.
- Antigen-quantity core is dead at full power (416 samples, perm p=0.77).

**Implication for any new modeling work:** do not re-open the representation
question on optimistic priors. Every new predictor must be benchmarked LOCO
against the immune-composition floor and reported honestly even if it loses.

## Active branch: open two-branch predictor vs GEM-1 (session `15defe54`)

**Goal:** build an *open* bulk-RNA ICB-response predictor to benchmark against
Synthesize Bio's (closed) GEM-1, using the raw reads. Plan approved 2026-07-08.

**Design (why it is not just "another learned representation"):**
- **Branch 1 — portable expression latent.** scVI / scGPT frozen embeddings
  *anchored on a large public bulk corpus* (recount3 SKCM + TCGA-SKCM + GTEx) so
  latents transfer across cohorts. This is the open analog of GEM-1's mechanism
  (batch-denoised portable latents) — aimed squarely at the transfer failure the
  kill-test found.
- **Branch 2 — non-reference raw-read block.** The signal GEM-1 structurally
  cannot see: RNA editing (AEI, done n=16), intron retention, splicing, and
  **TE/ERV (the open gap — Telescope parse-verified, never run on real data)**.
  Additive feature block on the leakage-guarded head.
- Honest head = `src/model.py::response_organization` (already fold-contained,
  LOCO + comparator + incremental AUROC). Rigor via `signature-rigour-harness`.
- Head-to-head on **IMvigor210** (open processed expr; a GEM-1 validation cohort)
  for a directly comparable AUROC.

**Autonomous decisions in effect (override if wrong):**
1. Raw-read expansion streams prefetch→salmon→delete-FASTQ per sample; hard floor
   ~150 GiB free (repo already 339 GB; disk has 360 GiB free). Balanced subset if
   a cohort would blow the budget.
2. GPU encoders (Evo2 / HyenaDNA / Caduceus fine-tune) stay **deferred** per
   `docs/ENCODER_REVIEW.md` — CPU-only 64 GB Mac. Branch 1 is CPU scVI/scGPT.
3. IMvigor210 is counts-only → tests Branch 1 only (fair vs GEM-1).
4. If the two-branch predictor does not beat the immune floor at LOCO, that is
   the reported result. No tuning-to-taste.

## Shared substrate produced by this branch (use it, don't rebuild)

`results/predictor/frozen_analysis_set.parquet` (+ `_meta.json`) —
149 labeled PRE samples (gide 73 / riaz 49 / hugo 27; R=65 / N=84), keyed on
`run_accession`, `loco_fold=cohort`. Carries: response, iAtlas proxy neoantigen
features, TMB/purity/ploidy/clonality, AEI editing (n=16), and a
`has_rawread_expr` flag (67 samples have de-novo salmon TPM today).
Label source = `run_catalog.resp_NR` (authoritative), PRE only.

## What other agents can do to help

- **TE/ERV (Telescope) on the staged HISAT2 BAMs** is the single most valuable
  open contribution — it is the branch-2 differentiator and has never run on real
  data. BAMs are in `results/editing_bams/` (16 samples, ~3 GB each).
- Expanding `results/features/quant_gene_tpm.parquet` beyond 67 raw-read samples
  directly grows branch 1's training set. Append on `run_accession`.
- When **liu2019** expression lands, the auto-runner (`analysis/kill_tests/
  move2_autorun.py`) fires the powered adjudication — coordinate so we re-run the
  two-branch LOCO at n~180 / 4 held-out cohorts too.

## Conventions

- Sample key is **`run_accession`** everywhere; join labels via the crosswalk
  rules in the `icb-rna-pilot-ingest` skill (never merge on [run, cohort]).
- Git metadata lives in `.gitmeta/` — drive with
  `GIT_DIR=.gitmeta GIT_WORK_TREE=. git ...`. `results/` is git-ignored.
- Feature-registry writes must use `analysis.registry_update.register_features`
  (atomic, flock) — concurrent sub-agents will clobber a plain `json.dump`.
- **Disk: use ARCHIVAL CRAM, not native-CRAM working format.** To reclaim disk
  on alignment BAMs, run `pipelines/scripts/archive_bam_to_cram.sh <ref.fa> <bam|dir>`
  AFTER the subworkflows have consumed the BAM: it does a lossless
  `samtools view -C -T <ref>`, deep-verifies the round-trip (TLEN-sign-normalized),
  and deletes the BAM. ~53% smaller, verified lossless on the pilot BAM. Do NOT
  convert the spine to CRAM-native — that would force the reference through
  featureCounts/TEcount and add encode/decode CPU on every pass for no extra
  saving (the BAM is transient anyway). See RUNBOOK §8.
- **TE/ERV (Telescope) now runs on real data end-to-end (pilot, exit 0).**
  `pipelines/te_erv/` runs bowtie2 multimap + Telescope EM (locus) + TEcount
  (family) natively on arm64. Telescope v1.0.3 is vendored as a patched arm64
  wheel (`pipelines/te_erv/vendor/`, rebuild via `build_telescope_wheel.sh`); run
  with `bash pipelines/env.sh` sourced and the te_erv samplesheet (fastqs + spine
  BAM). Pilot outputs published in `results/te_erv_integration_test/matrices/`:
  `telescope_locus_counts.tsv` (15,209 loci) and `tetranscripts_family_counts.tsv`
  (1,328 TE families + 63,140 genes); bowtie2 align 75.4% overall, EM converged
  (200 iters, log-likelihood ~2.04e7). See RUNBOOK §3 status + §5 items 14–15.

_Last updated: 2026-07-08 by session 15defe54 (Phase 1, frozen set complete)._
