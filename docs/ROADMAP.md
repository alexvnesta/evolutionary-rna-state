# Roadmap

Working plan; refined as we go. The falsifiable claims (see `README.md`) drive
the ordering: internal co-variation is tested before anything touches response
labels, and rigor guardrails are built in from the start rather than retrofit.

## Stage 0 — Organization (this stage)
- Fresh repo, compliance boundary documented, public-data metadata restored.
- Data inventory + run catalog understood and committed.

## Stage 1 — Data acquisition
- Pull open melanoma ICB raw reads (Gide 2019 `PRJEB23709`, Riaz 2017
  `GSE91061`/`PRJNA356761`) from ENA/SRA to git-ignored paths, verifying MD5
  against `data/manifests/download_receipts.json`.
- Reconcile run catalog ↔ clinical labels; freeze the pretreatment analysis set.

## Stage 2 — Interpretable feature backbone *(in progress)*
- Quantify transcript/gene abundance (decoy-aware pseudoalignment).
- Two quant arms in flight: an arm64 Nextflow pipeline (sibling session) and a
  lightweight de-novo salmon pilot (`docs/PILOT_NOTES.md`). The alignment step
  was **switched from STAR to HISAT2** because the only osx-arm64 STAR conda
  build (2.7.11b) produces empty BAMs on this toolchain. HISAT2 pilot verified:
  sample PD1_35_PRE aligns to a valid BAM with 79,080,188 reads, 100% mapped,
  100% properly paired (full-depth) — confirming the pivot resolves the
  empty-BAM failure.
- The leakage-guarded latent-state + response backbone (`model.py`) is already
  built and validated on mock data (latent recovery corr ≈ 0.96; response
  AUROC S = 0.833 vs TMB = 0.559).
- Score established immune/expression signatures.
- Estimate composition + purity (immune/stromal deconvolution) as **explicit
  covariates** the latent state must beat or be shown orthogonal to.

## Stage 3 — RNA-phenotype quantification *(underway)*
- Splicing, intron retention, RNA editing, TE activation, fusions, cryptic ORFs.
- Each phenotype quantified as independently as possible to avoid manufacturing
  co-variation (proxy-circularity guard).
- **RNA editing: real output on 16 samples** across both cohorts (Gide `ERR22089xx`,
  Riaz/Hugo `SRR*`) — 16 HISAT2 BAMs + 16 A-to-I editing-index (AEI) tables in
  `results/editing_bams/`. Values pass the signal/noise guard (e.g. ERR2208909
  AEI 0.226%, S/N 13.8 vs 0.016% noise floor).
- Intron-retention and splicing arms integrated onto the HISAT2 spine (code
  complete). TE/ERV (Telescope) parse-verified but not yet run on real data.

## Stage 4 — Latent state (internal claim)
- Learn a low-rank / unsupervised representation of the RNA phenotypes.
- Test co-variation **without** response labels.

## Stage 5 — Organization of response (external claim)
- Relate the latent representation to ICB response.
- Compare against / adjust for TMB, purity, composition. Split sits **above**
  any encoder or feature-selection step (CV-leakage guard).
- External validation on the held-out cohort.

## Stage 6 — Raw-read / sequence-model encoder branch (decision made; build deferred)
Complementary branch that embeds sequence directly so non-reference signal
contributes without committing to any one annotation. The "fine-tune vs.
off-the-shelf embed; which model" question has been researched and resolved —
see `docs/ENCODER_REVIEW.md`. Current decision, in build order:

- **Frozen embeddings first, not fine-tuning.** At n ≈ 150 melanoma
  pretreatment samples on CPU-only hardware, the two 2025 DNA-foundation-model
  benchmarks (Tang et al.; Feng et al.) predict that fine-tuning a raw-sequence
  genomic LM as the primary encoder would overfit and lose to the existing PCA
  backbone. Extract embeddings from a frozen HyenaDNA or Caduceus (state-space,
  single-nucleotide resolution, long context, RC-equivariant — the right fit;
  attention models' quadratic context is the wrong one) and feed them into the
  existing leakage-guarded latent-state and response models as an added feature
  block and comparator. The low-data regime is exactly where frozen embeddings
  are most likely to help.
- **Parameter-efficient fine-tuning (LoRA/adapters) is gated**, not cancelled.
  Undertake it only when (a) frozen embeddings show signal a linear head
  underuses AND (b) a GPU host is available (not present yet). When it happens,
  aim the training signal at **phenotype reconstruction** (mapping embeddings
  onto the six RNA-phenotype axes), never directly at ICB response — that
  separation is what the falsifiable claims require and avoids the small-n trap.
- **No encoder work is active or scheduled** until those gates are met. Nothing
  is currently fine-tuning or embedding reads; the raw-read work in flight is
  Stage 2/3 quantification (below), not an encoder.

## Rigor guardrails (applied throughout)
Composition confounding · proxy circularity · CV leakage · provenance.
The `signature-rigour-harness` skill is loaded at the modeling stages.

## Related design docs
- `docs/ENCODER_REVIEW.md` — encoder choice and the fine-tuning decision (Stage 6).
- `docs/PILOT_NOTES.md` — salmon de-novo quant pilot, arm64 toolchain notes.
