# Feasibility protocol: bulk-RNA inference of clonal-evolutionary trajectory class as an ICB-response predictor

## Purpose and design

Learn a discrete library of clonal-evolutionary **trajectory classes** from tumors with matched single-cell RNA-seq and DNA, then train a classifier that assigns a tumor to its class from **bulk RNA alone**, and test whether that class predicts immune-checkpoint-blockade (ICB) response. The trajectory class is a *supervised label* transferred into bulk RNA; this is label-recovery, not de-novo phylogeny reconstruction.

## What the data scan established (hard numbers)

No public deposit carries all four layers (scRNA + WXS + bulk RNA + ICB outcome) in the same patients, so the pipeline is federated across three tiers. Matched-patient counts were derived by parsing GEO sample sheets (patient ID + assay type per sample); where a sheet did not encode patient identity, the count is marked not-determinable rather than guessed:

| Tier-1 cohort | Accession | Matched patients | Layers | Class |
|---|---|---|---|---|
| cSCC | GSE144240 + GSE144237 | **9** (9/10 with both scRNA & WES) | scRNA + WES + spatial | solid |
| DDLPS | GSE221492 + GSE221493 | **16** (matched sc + bulk) | scRNA + bulkRNA (+WGS in paper, not in these GEO series) | solid |
| AML | GSE156934 | **38** (parsed patient IDs; scDNA study) | single-cell DNA | heme |
| CLL | GSE161610 (sc) + GSE161711 (bulk) | **5** (matched sc + bulk) | scRNA + bulkRNA | heme |
| mel/NSCLC brain | GSE192402 + GSE216069 | n.d. | scRNA + WGS | solid |

Counts marked with a bold number were computed from the GEO sample sheets by parsing patient IDs and assay type; **n.d.** (mel/NSCLC) could not be resolved because the sample titles are assay descriptors with no per-patient identifier. Two caveats corrected from an earlier draft: CLL's matched **scRNA+bulk** overlap is only **5** patients (the single-cell series GSE161610 covers 5 patients; the paper's "29" refers to patients matched across three *bulk* sample types within GSE161711, not sc+bulk), and AML's 38 is a parsed patient-ID count from a single-cell-**DNA** study — a per-patient DNA+RNA co-assay was not confirmed from the sheet.

Training budget: **~25 solid-tumor** matched patients (cSCC 9 + DDLPS 16), rising to roughly **60–68** if the AML single-cell-DNA cohort (38) is used for clone-tree construction. The CLL sc+bulk overlap (5) is too small to contribute meaningfully. Anchor tier (scRNA + ICB response): GSE120575, GSE115978, GSE123813, GSE121636. Test tier (bulk RNA + ICB outcome): GSE91061 (Riaz, n=109, with matched WES), GSE78220 (Hugo, n=28), GSE115821 (Auslander, n=37), plus Gide (ENA) and Van Allen (dbGaP).

## What the power analysis established (calibrated go/no-go)

A 3-class recovery simulation, with the noise floor calibrated to the real gene-gene covariance of the Hugo melanoma bulk FPKM matrix (GSE78220, 28 samples x 25,268 genes), gives balanced-accuracy as a function of training N and effect size (how strongly the clonal class imprints bulk expression; chance = 0.33, usable threshold ~0.70):

- At **N=25** (solid-tumor budget): 0.51 at moderate imprint, 0.73 only at strong imprint. **Underpowered unless the class signal is strong.**
- At **N≈65** (adding the AML clone-tree cohort): 0.82 moderate, 0.97 strong — **comfortably usable if a moderate signal exists.**
- At weak imprint (heavy bulk dilution): 0.53 at N=67 and 0.64 at N=100, crossing the ~0.70 usable line only around N=150 (0.79). So a weakly-imprinting class is recoverable, but only with a training set roughly double what the current public matched cohorts provide.

The decisive contingency: feasibility hinges on **effect size**, i.e. whether a clone-defined trajectory class leaves a recoverable imprint on bulk average expression. That is itself an empirical quantity, and Phase 0 below measures it before any classifier is built.

## Phase 0 — Measure the effect size first (1-2 weeks)

The cheapest decisive experiment. Do not build the full pipeline until this returns a number.

1. Take one Tier-1 cohort with clean clones (cSCC GSE144240/237, or AML GSE156934 for the cleanest trees).
2. Assign cells to clones: cardelino (WES tree + scRNA variants) for cSCC; native scDNA genotypes for AML. Reconstruct the clone phylogeny with PyClone/Canopy on the DNA.
3. Pseudobulk each tumor by summing clone-resolved expression back to a whole-tumor profile.
4. Compute the between-class variance ratio in pseudobulk space for candidate class definitions (branched vs linear; truncal- vs subclonal-driver). This *is* the effect size the power curve is parameterized on. Map it onto the eff=0.6/1.0/1.5 axis.

Decision gate: if pseudobulk effect >= moderate (eff ~1.0), proceed to the heme-inclusive N=67 design. If weak, restrict to a binary (2-class) label, or escalate to controlled-access cohorts (TRACERx, dbGaP) for larger N before proceeding.

## Phase 1 — Clone assignment and trajectory-class library (4-6 weeks)

1. Per Tier-1 cohort, run clone assignment (cardelino / Numbat where DNA is thin) and DNA-side subclone clustering (PyClone/ABSOLUTE).
2. Reduce each tumor's tree to a feature vector: branching index, truncal-driver identity, clonal vs subclonal neoantigen fraction (netMHCpan on WES-called mutations x cancer-cell fraction), HLA-LOH status (LOHHLA).
3. Cluster tumors into K classes (start K=3, TRACERx/Liao precedent; REVOLVER for cross-patient recurrent trajectories). Fix K by silhouette + biological interpretability, not by fit alone.
4. Deliverable: a labeled table (patient -> trajectory class) plus the class centroids in tree-feature space.

## Phase 2 — Bulk-RNA classifier (3-4 weeks)

1. Training features: pseudobulk the scRNA (data augmentation) + real matched bulk where it exists (DDLPS, CLL). Elastic-net or gradient-boosted classifier; nested CV; report balanced accuracy with the N-appropriate CI from the power grid.
2. Explicitly model the mixture: include tumor purity (ESTIMATE/ABSOLUTE) as a covariate, since bulk dilution is the dominant confound and the reason weak-imprint classes fail.
3. Benchmark honestly against the scalar-ITH baselines (SpliceHetero, transcriptomic-ITH) and a pathway-signature baseline (Stem.Sig): the novelty claim requires beating "just a heterogeneity score."

## Phase 3 — ICB-response test (2-3 weeks)

1. Apply the frozen classifier to the bulk ICB cohorts (Riaz, Hugo, Auslander; Gide/Van Allen on access).
2. Primary endpoint: association of predicted trajectory class with RECIST response and PFS/OS (multivariable, adjusting for TMB and purity).
3. Pre-register the hypothesis direction from the biology: truncal-/clonal-neoantigen-dominant classes should enrich for responders (McGranahan); subclonal-diversified / HLA-LOH classes for non-responders.

## Risks and mitigations

- **Small, heterogeneous training N.** Mitigate by pooling cancer types for the *class model* while keeping the bulk classifier cancer-aware; lean on heme for clean trees; hold TRACERx/dbGaP as the escalation path. The power grid says N=67 is the realistic floor for a usable 3-class model.
- **Solid-tumor phylogenies are noisier than heme.** Numbat's haplotype-aware CNV phylogeny is the fallback where SNV coverage in scRNA is too sparse for cardelino.
- **Cross-cancer label leakage.** A trajectory class must not be a proxy for tissue of origin; test by within-cancer-type cross-validation.
- **Purity confound.** Carried as an explicit covariate in Phase 2, not ignored.

## Bottom line

Every ingredient is public and reachable. The design is sound and squarely on the project thesis. The single quantity that decides success — how strongly a clonal-trajectory class imprints bulk expression — is measurable in Phase 0 for one to two weeks of work, and the calibrated power curve says the difference between "restrict to a binary label / go controlled-access" and "3-class model works at N=67" turns entirely on that number. Measure it first.
