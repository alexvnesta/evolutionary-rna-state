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

## Stage 2 — Interpretable feature backbone
- Quantify transcript/gene abundance (decoy-aware pseudoalignment).
- Score established immune/expression signatures.
- Estimate composition + purity (immune/stromal deconvolution) as **explicit
  covariates** the latent state must beat or be shown orthogonal to.

## Stage 3 — RNA-phenotype quantification
- Splicing, intron retention, RNA editing, TE activation, fusions, cryptic ORFs.
- Each phenotype quantified as independently as possible to avoid manufacturing
  co-variation (proxy-circularity guard).

## Stage 4 — Latent state (internal claim)
- Learn a low-rank / unsupervised representation of the RNA phenotypes.
- Test co-variation **without** response labels.

## Stage 5 — Organization of response (external claim)
- Relate the latent representation to ICB response.
- Compare against / adjust for TMB, purity, composition. Split sits **above**
  any encoder or feature-selection step (CV-leakage guard).
- External validation on the held-out cohort.

## Stage 6 — Raw-read encoder branch (design under discussion)
- Complementary branch embedding reads directly so non-reference signal
  contributes. Scope (fine-tune vs. off-the-shelf embed; which model) is an
  open design question to resolve before build.

## Rigor guardrails (applied throughout)
Composition confounding · proxy circularity · CV leakage · provenance.
The `signature-rigour-harness` skill is loaded at the modeling stages.
