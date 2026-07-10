# NORTH STAR — the one thing this project is finishing

_Set 2026-07-10 by session `23cf8106` (user-directed), for the Built-with-Claude hackathon (ends 2026-07-13)._

## The hypothesis (verbatim, user)
> The captured bulk RNA-sequenced transcriptome of a human tumor contains sufficient information
> to reconstruct its latent evolutionary state in the context of immunotherapy response.

The project explores a **hybrid** approach: learned representations from raw RNA-seq reads **+**
interpretable non-reference transcriptomic features (alternative splicing, intron retention, RNA
editing, TE/ERV activity, fusion transcripts) to identify signatures of ICB sensitivity.

## The single goal for the deadline
**Produce the first cohort-scale NON-REFERENCE feature matrix and run one pre-registered two-block
LOCO test: non-reference features vs the immune floor.** Everything else is subordinate to this.

Why: a four-audit forensic synthesis (`AUDIT_SYNTHESIS.md`, `756ba1cd`) established that every prior
"DEAD" verdict rested on gene expression, expression-PCA, or a DNA/WES antigen-quantity proxy — **never
on the non-reference transcriptome features the hypothesis is actually about.** The central claim is
**largely UNTESTED, not falsified.** The non-reference features are 0/106 built at cohort scale. This is
the gap; closing it is the north star.

## What is IN scope (finishable, novel vs COMPASS)
1. **Unified Apple-Silicon Nextflow workflow**: one top-level DSL2 workflow, HISAT2 spine → fan-out to the
   five validated per-feature subworkflows (rna_editing, te_erv, intron_retention, rnasplice, rnafusion)
   → one merged per-sample feature matrix. A wrapper of proven arm64 arms, not a rewrite.
2. **First cohort-scale non-reference feature matrix** keyed to the authoritative n=106 frame
   (`reconciled_frame_n106.parquet`). Target n=106; staged fallback = whatever is aligned.
3. **One pre-registered two-block LOCO test** (non-ref block vs immune floor), with the four
   signature-rigour checks baked in (proxy circularity, composition confounding, CV leakage, provenance).
4. **Dated hackathon result brief + figures + updated PROJECT_STATUS**, honest verdict either way.

## What is OUT of scope for the deadline (scoped as future work)
- **The learned-representation / encoder half of the hybrid.** The forensics showed no encoder was ever
  run on patient-specific/aberrant sequence: EVA collapsed to a fixed linear map of expression (R²=1.000),
  Orthrus never ran at patient scale. Doing it correctly requires building personalized per-patient
  transcript sequences (reference + somatic variants + retained introns + novel junctions + edited sites +
  TE loci + fusions), embedding those, and pooling within-patient — weeks + GPU. A rushed version would
  reproduce the `Wn @ E` expression-collapse that made the last run meaningless. Delivered as a **protocol**
  (`docs/ENCODER_PHASE_PROTOCOL.md`), NOT a claimed result.

## Honest status of the substrate (verified on disk 2026-07-10)
- Aligned HISAT2 BAMs: **32** unique in `results/rnaseq_cohort/hisat2/` (74 of 106 remain).
- **All 30 label-mapped BAMs are gide2019** (16R/14N). Precomputed non-ref caller outputs
  (`results/nonref_run/out/`, 8 samples) add 6 Gide + 2 Riaz. **0 Hugo aligned.**
- CONSEQUENCE: the immediately-computable result is a **within-Gide cross-validated** non-ref-vs-floor
  test. A genuine **LOCO** (leave-one-cohort-out) needs ≥2 cohorts at usable n — so the background
  alignment MUST prioritise **Hugo (27 samples, cheapest ~239 GB)** and more **Riaz** to unlock LOCO.

## Compute path (accurate)
- **Modal (`byoc:modal`) is available but has NO env image built yet** — a first-time image build is
  required (via the `compute_provider` kernel / `compute-env-setup` skill) before cloud fan-out can run.
  This is a not-yet-built state with a documented build path, NOT a broken provider.
- Local batched alignment via `pipelines/scripts/run_cohort_batched.sh` with the validated AWS
  fast-fetch (~3.5 min/sample fetch) is the fallback and works today.

## Coordination
- Sole coordinator: `837512d2` (owns `PROJECT_STATUS.md`).
- This finalization is driven by session `23cf8106` directly against the shared repo (sibling root
  sessions cannot be messaged across the comms topology; steering is via these docs + direct ownership).
- Encoder/hybrid work is future-work; do not present it as a deadline result.
