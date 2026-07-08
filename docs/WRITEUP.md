# Is there a latent evolutionary RNA-state that organizes immune-checkpoint response?

**A rigorous test on five public melanoma ICB cohorts — and where the signal actually lives.**

*Built with Claude: Life Sciences hackathon, Research Track. All code authored during the event; all data public.*

---

## Abstract

Response to immune-checkpoint blockade (ICB) is only weakly predicted by the
approved biomarkers. On five harmonized public melanoma ICB cohorts we
reproduce that ceiling: tumor mutational burden (TMB), the TIDE signature, and
a clinical panel all sit at pooled cross-validated AUROC 0.53–0.62. The
project's hypothesis is that early driver mutations set an evolutionary
trajectory whose downstream **RNA phenotypes** — alternative splicing, intron
retention, RNA editing, transposable-element (TE) activation, fusion
transcripts — are coordinated manifestations of a single latent *RNA-state S*
that shapes antigenicity and response. We test two falsifiable predictions.
**(1) Internal:** the RNA phenotypes should share low-rank structure even
without response labels. Using the WES-derived neoantigen proxies available at
scale (iAtlas), this is **null** (permutation p = 0.78; the only strong
co-variation is between two burden-derived proxies, which a circularity guard
attributes to shared mutational burden). **(2) External:** a representation of
S should stratify response beyond TMB. The WES-proxy phenotypes add
**zero** incremental AUROC over the TMB/TIDE floor, and a power analysis shows
the proxy carries an implied effect ~10× too small to detect at this n.

We argue the null is diagnostic, not disconfirming: WES-derived proxies cannot
see a transcriptomic state. We therefore built a lightweight raw-read pipeline
from scratch (salmon, arm64-native) and de-novo quantified the transcriptome of
**40 pre-treatment samples** (Gide 2019 n = 30; Riaz 2017 n = 10) directly from
ENA FASTQ. A de-novo antigen-presentation axis (B2M/HLA-A/TAP1) separates
responders **within** Gide (leave-one-out AUROC ≈ 0.87), but rigorous analysis
shows the signal is **largely infiltration-driven** (tumor-intrinsic residual
AUROC ≈ 0.62) and **does not replicate** on held-out Riaz (AUROC 0.36–0.44).
The honest conclusion: a single-phenotype, gene-level de-novo axis reflects the
microenvironment and is cohort-specific — *not* a transferable tumor-intrinsic
RNA state. Testing the actual hypothesis requires the multi-phenotype,
WES-blind features (intron retention, RNA editing, TE/ERV) that standard
pipelines discard; we built and staged the infrastructure to produce them.

---

## 1. Background and hypothesis

Every clinically proven ICB biomarker is a DNA-level feature (dMMR/MSI-H,
TMB-high), a bulk-expression signature (IFN-γ / T-cell-inflamed GEP), a
single-protein stain (PD-L1 IHC), or an external factor. None captures the
**non-canonical RNA state** of the tumor: splicing-derived neoantigens,
TE-derived antigens, intron retention, RNA editing, and fusion transcripts.

The project hypothesis is that these are not independent biomarkers but
observable manifestations of one underlying **evolutionary RNA-state S**. Early
driver mutations establish a trajectory; over time genomic instability,
epigenetic remodeling, RNA-processing dysregulation, and immune selection
produce *coordinated* transcriptomic abnormalities that determine tumor
antigenicity and, ultimately, response to ICB.

This yields two falsifiable claims:

- **Internal (no labels):** the RNA phenotypes co-vary — they share low-rank
  structure — beyond what mutational burden alone explains.
- **External (with labels):** a representation of S stratifies ICB response
  *beyond* TMB, tumor purity, and immune composition.

## 2. Data and the baseline ceiling

Five public melanoma ICB cohorts, cBioPortal/iAtlas-harmonized: Gide 2019
(PRJEB23709), Riaz 2017 (GSE91061), Hugo 2016 (GSE78220), Liu 2019, and DFCI
2019. After freezing to one pre-treatment baseline per patient, the analysis
set is **416 samples**; the phenotype-bearing subset (cohorts with iAtlas
neoantigen categories) is **264**.

The approved-biomarker ceiling, reproduced with pooled stratified-group
cross-validation:

| Model | Pooled CV AUROC | 95% CI | n |
|---|---|---|---|
| TMB (nonsynonymous) | 0.62 | 0.54–0.70 | 190 |
| TIDE signature | 0.53 | 0.47–0.60 | 270 |
| Clinical panel (TMB+MUT+ICI) | 0.62 | 0.54–0.70 | 190 |

This weak floor is the gap the RNA-state hypothesis proposes to fill. Any
RNA-derived feature must beat it to be claimed as adding signal.

## 3. The internal claim is null on WES-derived proxies

The neoantigen categories available at scale (iAtlas SPLICE / ERV / FUSION)
are **WES-derived proxies** — counts of predicted neoepitopes from exome
variant calls, not direct transcriptome measurements. Testing co-variation on
these:

- Pairwise partial-Spearman between splice and fusion proxies, residualized on
  mutational burden and cohort: pooled ρ = **−0.02** (p = 0.80); per-cohort
  signs are inconsistent.
- Permutation low-rank test on the burden-residualized splice+fusion block:
  observed PC1 variance fraction **0.510** vs null mean 0.529 (95th percentile
  0.571), **permutation p = 0.78** (n = 192, 5,000 permutations). No shared
  low-rank structure.
- The one strong off-diagonal in the full phenotype matrix is SNV–INDEL
  (ρ = 0.50). A circularity guard shows SNV-neoantigen shares R² = 0.53 of its
  variance with burden (|Spearman| = 0.85) — i.e. it *is* a burden proxy, which
  is exactly the trivial co-variation the test is designed to strip.

An independent rigor harness (four guards: circularity, cohort confounding,
CV-leakage, missingness) corroborates: the SPLICE proxy's naive
response-association (OR 1.34, p = 0.03) **does not survive cohort adjustment**
(OR 1.13, p = 0.42) — a cohort artifact, not biology; ERV is degenerate
(nonzero in < 2% of samples).

## 4. The external claim adds nothing — and the power is too low anyway

Over a TMB floor of pooled-CV AUROC 0.61 (n = 190), adding each WES-proxy
phenotype yields **incremental AUROC indistinguishable from zero**: splice
ΔAUROC −0.020 [−0.094, +0.055], fusion −0.019 [−0.079, +0.034], TE/ERV −0.011
[−0.038, +0.010], indel −0.016 [−0.069, +0.038], SNV +0.018 [−0.013, +0.052]
(point estimate + paired-bootstrap 95% CI; every interval straddles zero). A
meta power analysis makes the null interpretable rather than
merely underpowered-ambiguous: at n = 192 this design can detect a variance
explained of **≥ 22.6%** at 80% power (≥ 15.3% at 50%), but the observed
low-rank structure implies only **~2.0%** — roughly an order of magnitude below
the detection floor. The proxies do not carry the effect, and we could not have
seen it if they did.

**Interpretation.** This is the pivot of the project. A null on a WES-derived
proxy is not evidence against a *transcriptomic* state — it is evidence that
exome-based annotations cannot test the hypothesis. Splicing, intron retention,
RNA editing and TE activation are RNA-level phenomena; the variant-calling and
annotation pipelines that produce iAtlas proxies discard exactly the
non-reference signal where the fingerprint would live. To test the hypothesis
one has to go back to the **raw reads**.

## 5. A raw-read arm: de-novo transcriptome quantification from scratch

With no GPU and a ~1.4 MB/s ENA bandwidth wall, a full 1.3 TB FASTQ pipeline
was infeasible in the hackathon window. We built a bandwidth-aware pilot:
salmon 2.3.1 (Rust/piscem, arm64-native), a GENCODE v44 transcriptome index
(251,955 transcripts), and a stream-subsample / stream-align-delete loop that
takes the first 3M read-pairs of each sample, quantifies, and deletes the FASTQ
— keeping peak disk bounded and finishing in minutes per sample. Every download
is gzip-integrity-checked before quantification; the guard correctly rejected
transiently-truncated downloads rather than feeding salmon corrupt data.

We de-novo quantified **40 pre-treatment samples** (Gide n = 30, 16R/14N; Riaz
n = 10, 5R/5N) at 81–95% mapping rates (range 80.8–95.4%), collapsing 251,955 transcripts to
62,266 genes per sample. From these purely de-novo features we built curated
immune/antigenicity axes (antigen presentation: B2M/HLA-A/TAP1/TAP2/HLA-B/PSMB9;
IFN mimicry; cytolytic) and a pan-leukocyte **infiltration** score
(PTPRC/CD3/CD8A/…) as the confound to control.

## 6. Within-cohort signal — but infiltration-driven and non-transferable

All cross-validation is **fold-contained**: panel z-standardization and any
covariate residualization are re-estimated on training samples only within each
fold, so no test-sample information touches feature construction.

1. **Within Gide (n = 30):** the antigen-presentation axis separates response at
   fold-contained LOO AUROC ≈ **0.87**, modestly above an infiltration score
   (0.80) and IFN mimicry (0.75).
2. **Tumor-intrinsic test:** after residualizing the immune-infiltration score
   out of the axis within each fold, the residual falls to LOO ≈ **0.62–0.64**.
   Most of the response signal tracks *how much immune infiltrate is present*,
   with only a weak tumor-intrinsic component. This is not an assertion from a
   single residualized number: the antigen-presentation axis is Spearman
   ρ = **0.77** correlated with the infiltration axis (R² = 0.57 shared
   variance), and *every* constituent antigen gene is individually
   infiltration-correlated — HLA-B 0.79, TAP1 0.78, B2M 0.71, HLA-A 0.65,
   PSMB9 0.62, TAP2 0.54. The de-novo "antigenicity" axis is largely reporting
   leukocyte content, which is why removing it collapses the response signal.
3. **Held-out replication (the decisive test):** an axis+logistic model trained
   on Gide and applied to held-out Riaz scores AUROC **0.36–0.44** — worse than
   chance. Within Gide the axis separates cleanly (R median +0.70 vs N −0.59,
   AUROC 0.91); within Riaz the responder and non-responder distributions
   overlap (R +0.32 vs N +0.11, AUROC 0.40, indistinguishable from chance at
   n = 5 vs 5). Standardizing within Riaz using its own statistics does not
   rescue it, so this is genuine non-replication, not a batch/scaling artifact.
4. **A second held-out cohort confirms it is not Riaz-specific.** We de-novo
   quantified 12 additional Hugo 2016 pre-treatment samples (6R/6N, 86–98%
   mapping) as an independent second held-out set. The Gide-trained axis scores
   AUROC **0.58 on Hugo** (95% CI [0.22, 0.92], n = 12) — better than on Riaz
   but still statistically indistinguishable from chance, and far below the
   within-Gide 0.87. Across **both** independent cohorts the axis fails to
   transfer, so non-replication is a general property of the single-phenotype
   gene-level axis, not a quirk of one cohort. The infiltration score is the
   only axis with any cross-cohort consistency, and even that is weak
   (Riaz 0.68, Hugo 0.47).

The n = 12 balanced pilot that opened this arm showed a perfect separator
(LOO AUROC 1.00); scaling to n = 40 with held-out validation revealed that as
small-sample optimism. Reporting the honest ≈ 0.87 within-cohort / non-replicating
number is the point of the deepening.

## 7. Rigor and robustness

Both arms were stress-tested rather than taken at face value.

**WES-proxy arm (four-guard harness).** A circularity guard flags SNV-neoantigen
as a mutational-burden proxy (R² = 0.53 shared variance, |Spearman| = 0.85) —
correctly identifying the one trivial co-variation in the matrix. A
cohort-confounding guard shows the SPLICE proxy's naive response association
(OR 1.34, p = 0.03) does not survive cohort adjustment (OR 1.13, p = 0.42): a
batch artifact, not biology. A CV-leakage guard on the clean pipeline collapses
to chance under label permutation (permuted AUROC ≈ 0.50), confirming no
information leak. ERV is degenerate (nonzero < 2%).

**De-novo arm (three robustness checks).** The within-Gide antigen signal is
(i) **leakage-free** — under fold-contained LOO with permuted labels the AUROC
collapses to a null mean of 0.29 while the observed 0.89 gives permutation
p = 0.005; (ii) **CV-scheme-stable** — LOO (0.888) and repeated stratified
5-fold (0.895 ± 0.009) agree within 0.01; (iii) **panel-definition-stable** —
every individual antigen gene gives 0.79–0.88 on its own, so the result is not
driven by one gene or one panel choice. The signal is genuine and robust
*within cohort* — which is precisely what makes its failure to transfer to two
held-out cohorts the substantive, defensible finding rather than a modelling
error.

## 8. Conclusion and forward path

On the evidence assembled here:

- The WES-derived neoantigen proxies show **no** low-rank RNA structure and
  **no** organized response signal beyond mutational burden — but they are the
  wrong instrument for a transcriptomic hypothesis.
- A de-novo, gene-expression-level antigen-presentation axis *does* separate
  response within a cohort, but it is **largely a microenvironment
  (infiltration) signal and does not generalize** across cohorts. Gene-level
  expression alone — essentially what standard bulk RNA-seq delivers — is
  insufficient to demonstrate a transferable tumor-intrinsic RNA state.

This sharpens rather than refutes the hypothesis. The actual test requires the
**WES-blind, multi-phenotype features** the standard pipelines throw away:
per-sample intron-retention ratios, RNA-editing (Alu editing index), TE/ERV
locus expression, per-sample splice-junction PSI, and fusion burden — quantified
from raw reads and assembled into the latent state S. We built and validated the
modeling machinery for exactly this (fold-contained latent-state construction,
permutation-null low-rank testing, leave-one-cohort-out loading stability, and
incremental-AUROC response organization over the TMB/TIDE floor), and staged the
raw-read infrastructure and a per-sample feature contract for the multi-phenotype
matrices. When those features land, S can be tested directly, with the honest
cross-cohort bar this project has already established.

## Figure guide (main-text arc)

The main figures follow the argument, not the chronology of the work
(8 figures, in `figure_deck.pdf`):

1. **Fig 1 — the punchline as hook.** The de-novo antigen axis separates
   responders within Gide (AUROC 0.91) and collapses to chance on held-out Riaz
   (AUROC 0.40): a candidate biomarker that looks real until tested honestly.
2. **Fig 2 — mechanism (why the standard approach fails).** WES-derived
   neoantigen proxies share no low-rank structure beyond mutational burden
   (permutation p = 0.78; co-variation matrices + permutation-null histogram).
3. **Fig 3 — the null is well-powered.** The co-variation test could detect a
   variance-explained of ~23% at 80% power, but the observed structure implies
   only ~2% — so the proxies are the wrong instrument, not merely underpowered.
4. **Fig 4 — the bar.** Approved biomarkers (TMB, TIDE, clinical panel)
   reproduce the modest known ceiling (~0.62 AUROC).
5. **Fig 5 — the external null.** No WES-proxy phenotype moves the needle over
   that floor (all ΔAUROC CIs straddle zero); RNA-derived features require raw
   reads.
6. **Fig 6 — the promising candidate.** Moving to raw reads, a de-novo
   antigen-presentation program cleanly separates responders in a balanced
   12-sample pilot — the hope the paper then stress-tests.
7. **Fig 7 — why it fails (mechanism of the confound).** The antigen axis is
   ρ = 0.77 correlated with leukocyte content (R² = 0.57), and every antigen
   gene is 0.54–0.79 infiltration-correlated — the axis largely reports
   infiltration, not tumor-intrinsic antigenicity.
8. **Fig 8 — the teardown.** The Gide-trained antigen axis fails to transfer to
   *both* independent held-out cohorts (within-Gide 0.87 → Riaz 0.36 → Hugo
   0.58, CI includes chance): non-replication is a general property, not a
   single-cohort quirk.

The field-coverage/missingness heatmap (Supp Fig S1, `figure_deck_supplement.pdf`)
is reference material for the supplement, not the narrative.

## Methods (brief)

Cohorts harmonized via cBioPortal/iAtlas; one pre-treatment baseline per
patient. Baselines: pooled stratified-group CV logistic models with per-fold
preprocessing (leakage-guarded). Internal claim: partial-Spearman with
burden+cohort residualization and a 5,000-permutation PC1-variance-fraction
null. External claim: incremental AUROC over TMB/TIDE with paired-bootstrap CIs;
meta power via variance-explained simulation. Raw-read arm: salmon 2.3.1,
GENCODE v44, 3M read-pair stream-subsample, gzip-integrity-verified;
gene-level TPM via tx2gene collapse. All response-organization CV is
fold-contained. Code and per-analysis artifacts are in the repository.

*Reproducibility: `src/` and `analysis/` contain the pipeline and modeling
code; `results/` contains every statistic quoted here as a saved artifact.*
