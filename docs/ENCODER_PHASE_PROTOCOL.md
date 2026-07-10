# Encoder-phase protocol — the learned-representation half (scoped future work)

**Status: NOT a deadline deliverable.** Scoped here so the hybrid half is a designed experiment, not an
afterthought. The interpretable non-reference half was tested this session (`HACKATHON_BRIEF.md`, negative);
this document specifies what the *learned-representation* half must do to be a real test rather than a repeat
of the failures the forensic audit already caught.

**v2 (2026-07-10):** revised after an expert-panel review (ML-rigour, genomic-FM, tumour-immunology lenses).
The panel identified design flaws that would have made the first version either uninterpretable or a wasted
GPU run. Their corrections are integrated below and flagged `[panel]`. The architecture is now the
**post-hoc attribution** design the user chose (two branches + attribution read off a trained head); the
jointly-trained COMPASS-style concept-bottleneck version is explicitly **post-hackathon future work**.

---

## 0. Hypothesis of record (sharpened)

Raw tumour RNA-seq contains **implicit, not-yet-named** features — carried in *non-reference sequence* the
annotation discards — that predict ICB response **beyond** what reference-expression / immune-composition
models (COMPASS-class) already capture. The deliverable model takes RNA-seq and emits (a) an ICB-response
prediction, (b) **manual** feature weights (interpretable caller + immune-floor coefficients), and (c)
**implicit** feature weights (attribution over learned sequence features). This is *not* a COMPASS repeat:
COMPASS is a concept bottleneck over **reference expression** that renounced implicit discovery for
interpretability; this design keeps the interpretable branch **and adds a learned sequence branch**, then
measures how much the learned branch adds over the named one.

## 0a. Inferential regime — DESCRIPTIVE, not confirmatory `[panel: ML-rigour, decisive]`

With **3 cohorts** the cohort-level permutation null has only 3! = 6 arrangements, so the finest attainable
p-value is **1/6 ≈ 0.167 — p<0.05 is unreachable for the primary endpoint regardless of signal strength.**
Therefore:
- The primary analysis is **descriptive / hypothesis-generating**: report incremental ΔAUROC point
  estimates, Hanley-McNeil + bootstrap CIs, effect sign, and per-fold consistency. **No "significant"
  language.**
- Confirmatory significance testing is deferred to a **pre-registered replication once ≥4–5 cohorts exist**.
- This limitation is stated **up front**, not discovered at analysis time.

## Why it is deliberately deferred
The forensic audit (`ENCODER_EVALUATION_FORENSICS.md`) established two facts that make a rushed encoder run
worse than none:
1. **EVA collapsed to expression.** The per-patient EVA feature was a fixed linear map of gene expression:
   `EVA_SAMP = Wn @ E`, R² = 1.000000 reconstructing 1024 EVA dims from 39 expression PCs, PC1 r=0.957.
   A "learned representation" that is an invertible linear function of expression tests nothing new.
2. **Orthrus never ran at patient scale.** The only Orthrus execution was a 5-gene CPU benchmark
   (B2M/GAPDH/ACTB/CD8A/TP53). The cited "~0.49 Orthrus LOCO" was a misattributed PCA/scVI baseline. Zero
   encoders were ever run on patient-specific or aberrant sequence.

The single lesson: **a learned representation is only worth testing if it ingests information expression
does not contain** — i.e. actual (patient-specific, non-reference) *sequence*, not a re-encoding of the
reference expression vector.

## The one hypothesis this phase tests
> A sequence model that reads a tumour's *own* non-reference RNA sequence produces features that predict ICB
> response **beyond** the immune floor AND beyond the interpretable non-reference block AND beyond TMB/purity.

If the representation cannot beat the floor, the hypothesis is unsupported at the achievable n. If it cannot
beat floor+TMB+purity **in inflammation-matched tumours**, it is re-deriving inflammation, not discovering
novel antigenicity `[panel: immunology]`.

## 1. What may be fed to the sequence model — the SEQUENCE-VISIBLE partition `[panel: genomic-FM, decisive]`

A sequence-only model (Evo 2 / HyenaDNA) sees **nucleotides, not abundance**. Feeding it the *reference*
locus of an abundance/processing event gives **every tumour the identical sequence → a zero-variance feature**
dressed up as a predictor. Partition the event classes:

- **SEQUENCE-VISIBLE (feedable — the tumour actually carries non-reference bases):**
  RNA-editing sites (A→I read as A→G), fusion/novel-junction breakpoint contigs, somatic SNV/indel/SV
  neo-sequence. **The Evo 2 branch is built ONLY from these.**
- **NOT sequence-visible (Evo 2 adds nothing — leave to the manual/abundance branch):**
  intron-retention fraction, TE/ERV family CPM, cryptic-ORF *expression level*. These are abundance/processing
  changes; scoring their reference locus yields constants.

Mandatory: **log the cross-tumour variance of every candidate Evo 2 feature and hard-drop zero/near-zero
variance features before the head.** For any genomic-locus feature, assert the input sequence differs across
tumours or the feature is inadmissible.

## 2. Feature construction (named scalars, not high-dim embeddings) `[panel: ML + genomic-FM]`

To keep "implicit feature weights" both *estimable at n≈106* and *interpretable*, reduce Evo 2 output to
**per-event / per-event-class named scalars**, not a raw latent vector:

- **Per event:** background-normalised **delta-likelihood / surprisal** — score the tumour sequence AND a
  matched-normal/canonical background (GTEx skin + REDIportal for editing; GENCODE canonical junctions for
  splicing; parental native windows for fusions) through the **same model + window**, and use the tumour
  score's **percentile/z within that background**. This converts an out-of-distribution likelihood (which
  otherwise just measures "is this a junction") into an interpretable aberrant-vs-normal score.
- **Windows:** site/junction-centred DNA windows (~1–4 kb, within pretrain context) with the tumour's
  non-reference base(s)/junction substituted in — **not** the full mature mRNA (exon-exon junctions are OOD
  for a genome model).
- **delta-likelihood construction per class must be written down:** editing/SNV = clean autoregressive LLR
  over a site-centred window; fusion/novel junction = window mean log-likelihood vs parental + canonical
  background, labelled a *naturalness* score (no clean reference allele exists); reference construction fixes
  the sign, so it is pre-specified.
- **Per sample:** a compact fixed vector per class — `{depth-normalised event count, mean LLR, max LLR,
  abundance-weighted-sum LLR}`. That is the entire Evo 2 block: a handful of interpretable features.
- **Dose (abundance) enters here, saturating not raw** `[panel]`: use `Σ_e LLR_e · log1p(CPM_e)` (antigen
  presentation saturates with dose), and ALSO include the two main effects (aberrancy alone, abundance alone)
  so the head can reveal the interaction is unnecessary. Rank/quantile-normalise each feature across samples
  to kill cohort scale effects. Pre-specify the dose form (product/log/rank) or charge the choice to the null.

## 3. The presentation layer — the biology bridge `[panel: immunology, decisive]`

Evo 2 surprise measures **genomic naturalness/constraint**, which is closer to *deleteriousness* than to the
property ICB cares about: **foreignness to the T-cell repertoire**, HLA-restricted and downstream of
translation. Do NOT treat raw delta-likelihood as immunogenicity. Insert the missing middle steps:

1. **Translation/NMD filter:** require an in-frame ORF across the aberrant junction/region; apply NMD-escape
   rules (PTC >50–55 nt upstream of the last exon-exon junction → degraded). Retained introns and fusions are
   heavily NMD-regulated; unfiltered they mostly make no protein.
2. **Presentation:** call patient **HLA-I** from RNA-seq (arcasHLA/OptiType), run candidate 8–11mers through
   **NetMHCpan-4.1 / MHCflurry 2.0 for the patient's own alleles** → a *presented aberrant-peptide load*.
   This CPU/arm64-friendly layer is the physics-based bridge that makes the implicit attribution mean
   immunology; **build and validate it FIRST — it de-risks the chain independent of whether Evo 2 runs.**
3. Keep raw Evo 2 delta-likelihood only as a **secondary "aberrancy prior," explicitly not labelled
   immunogenicity.** Optionally add a peptide-space **self-dissimilarity / foreignness** metric
   (Łuksza-type) as a more faithful foreignness feature than genome-LM surprise.

## 4. Anti-collapse — REPLACE the in-sample R²>0.95 gate `[panel: ML + immunology, decisive]`

The original `R²(embedding ~ expression_PCs) > 0.95` gate is **broken**: computed in-sample at p≈n it hits
~1.0 under a true null (it measures degrees of freedom, not collapse). Replace with:
- **Primary gate = incremental ΔAUROC of the Evo 2-residual** (embedding orthogonalised to expression PCs,
  orthogonalisation fit **in-fold**) over the floor, LOCO, CI excluding 0. Collapse is only fatal if the
  residual carries no independent ICB signal — test the residual directly.
- **Held-out (cross-validated) reconstruction R²** kept ONLY as a diagnostic, reported per-cohort and pooled —
  never in-sample.
- **Sham-embedding negative control:** a random projection of expression to the same reduced dimension as the
  Evo 2 features. If the sham matches Evo 2's ΔAUROC, the Evo 2 signal is expression in disguise — a sharper
  collapse test than any reconstruction R².
- **Extend collapse tests beyond expression** `[panel: immunology]`: a sequence-varying feature can still
  collapse onto **TMB, tumour purity, library depth, IFNG/GEP**. Gate each Evo 2 feature (and the block
  jointly, via CCA/PLS R²) against `{expression PCs, TMB, purity, depth, IFNG}`.

## 5. Model, attribution, and stability `[panel: ML-rigour]`

- **Head = L2-penalised logistic regression**, penalty by nested in-fold CV. Standardised coefficients ARE
  the weights.
- **Do NOT run SHAP/IG/attention on a linear head** — for a linear model attribution collapses to
  standardised-coefficient × feature and adds no information; it is a garden-of-forking-paths multiplier with
  zero content. Attribution through a *frozen* Evo 2 attributes to the frozen model, not to anything trained
  on labels — skip it for this deliverable.
- **Report weight STABILITY as a first-class output:** bootstrap CI per weight, per-fold sign-consistency,
  cross-fold rank correlation of the weight vector. A weight whose sign flips across folds is reported as
  **unstable, not a finding.** (At n≈106 with LOCO, instability — not attribution value — is the real problem.)
- **Complexity budget:** total free parameters (features + pre-specified interactions) ≤ ~n/10–n/15 of the
  effective per-fold training n (~80–100). **Pre-specify interactions; do not search them.**
- **Everything in-fold:** PCA, scaling, event selection, orthogonalisation, dose-form choice — nothing that
  touches the held-out cohort may see it. Event *calling/selection* must use external annotation or in-fold
  data only (else double-dipping).

## 6. Covariates, controls, and label harmonisation `[panel: immunology + ML]`

- **Competing covariates the non-ref branch must beat, nested:** floor → floor+TMB → floor+TMB+purity →
  +non-ref, each with Hanley-McNeil CIs. Add **HLA-I heterozygosity + somatic HLA-I LOH** (both established
  ICB predictors — Chowell et al. 2018, *Science* 359:582, doi:10.1126/science.aao4572, which reports maximal
  HLA-I heterozygosity improving post-ICB survival AND HLA-I LOH associated with poor outcome) as
  covariates/stratifiers.
- **Harmonise the response label** across Gide/Hugo/Riaz (RECIST vs irRECIST; CR/PR vs SD/PD; durable
  benefit) and **model regimen/prior therapy** (Riaz includes anti-CTLA4-progressed; Gide mixes mono and
  ipi+nivo). Confirm pre-treatment biopsies only; record biopsy site.
- **Biological positive controls the pipeline MUST pass before any non-ref result is interpreted:**
  (1) recover TMB→response in melanoma; (2) recover the IFNG/T-cell-inflamed floor; (3) recover cytolytic
  activity (GZMA/PRF1); (4) recover an ERV/viral-mimicry score as *associated with the floor* (confirms the
  confound is measurable); (5) a technical control that Evo 2 delta-likelihood separates **ClinVar
  pathogenic vs benign** before trusting it on non-ref events.
- **Further negative controls, each through the identical pipeline:** within-cohort permuted labels (null);
  Evo 2 scores on reference/random matched loci (null; positive ⇒ leakage/batch); scrambled event→score
  pairing (signal killed); regress cohort+depth out of each Evo 2 score and report residual predictiveness.

## 7. Compute — Apple Silicon is NOT a path for Evo 2 `[panel: genomic-FM, decisive]`

Evo 2's stock implementation is **CUDA/FlashAttention-only; there is no MPS/Metal path** — any hour spent
forcing it onto the M-series workstation is wasted. Split the work:
- **Locally (arm64, MPS/CPU):** develop and debug the ENTIRE pipeline on **HyenaDNA** (autoregressive,
  hg38-trained, small) — validate the variance filter, background normalisation, and anti-collapse gates.
  Build the **HLA + NetMHCpan presentation layer** here too (CPU-friendly).
- **Final scoring pass:** **Evo 2-7B (bf16) on a single NVIDIA GPU via Modal** (`byoc:modal`, first-time env
  image build required — documented, not yet built). Escalate to 20B/40B (FP8, Hopper-class) only as a
  sensitivity check. The NVIDIA BioNeMo NIM hosted Evo 2 endpoint is a zero-setup alternative for
  likelihood/embedding scoring (few-kb input cap — verify current docs; compatible with the recommended
  windows). Inference-only forward passes — batch the windows, no training cost. Keep Evo 2 **frozen and
  version-pinned**; log the exact event set and windows for reproducibility and batch-confound audits.

## 8. Honest expectation `[panel: all three converge]`
The interpretable non-reference features carried no signal at n=25 (this session) and antigen *quantity*
carried none at n=106. Swapping in a frozen sequence model **adds parameters, not samples** — the a-priori
expectation is another (this time possibly *false*) negative, consistent with 2025 frozen-DNA-LM benchmarks.
Value is as a **correctly-designed, falsification-first** test: success is narrowly defined as **ΔAUROC CI
excluding 0 out-of-cohort on the sequence-visible classes, over floor+TMB+purity, in inflammation-matched
tumours** — with the outcome reported as effect size, not significance (§0a). Report which event class
(edits vs fusions vs junctions) drives any lift.

## 9. Concrete next steps (when resumed)
1. **Build the presentation layer first** (arm64, no GPU): arcasHLA on the cohort BAMs → NetMHCpan-4.1 on
   translated, NMD-escaping ORFs from the sequence-visible events → presented-peptide-load features.
2. Build and debug the full pipeline locally on **HyenaDNA** (variance filter, background normalisation,
   in-fold orthogonalisation, sham-embedding control, weight-stability bootstrap).
3. Assemble covariates: TMB, purity (ESTIMATE/ABSOLUTE), HLA het/LOH, harmonised response label, regimen.
4. Run the biological positive controls (§6) and the ClinVar technical control — **gate: all must pass.**
5. **Only then** build the Modal GPU image and run the frozen **Evo 2-7B** scoring pass on the
   sequence-visible events; add as block 3 to the (descriptive) LOCO harness with the nested-covariate ladder.
6. Report ΔAUROC + CIs + weight-stability; **no significance claims** until ≥4–5 cohorts (§0a).
