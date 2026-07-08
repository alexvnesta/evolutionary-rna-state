# Methods

## Design

### Overview and design

We test whether a tumor's clonal-evolutionary state imprints its transcriptome strongly
enough to be recovered from expression data by a supervised classifier, and whether that
state, read from bulk RNA, predicts immune-checkpoint-blockade (ICB) response. The design is
supervised label recovery: single-cell data is used to define clone-level transcriptional
states and a per-clone immune-state label, and the question is whether that label is
recoverable from pseudobulk (and, ultimately, real bulk). We do not attempt de-novo clonal
phylogeny reconstruction from bulk, which is underdetermined because bulk is an unknown
mixture of clones.

### Data sources

All data are public and were obtained from Gene Expression Omnibus (GEO) unless noted.
Single-cell cohorts, one per tumor lineage: cutaneous squamous cell carcinoma (cSCC,
GSE144236, with matched WES in GSE144237); dedifferentiated liposarcoma (DDLPS, single-cell
GSE221493 and matched bulk GSE221492); acute myeloid leukemia (AML, van Galen GSE116256);
and lung adenocarcinoma (LUAD, Kim GSE131907). Bulk ICB cohorts (melanoma): Riaz nivolumab
(GSE91061, with matched WES), Hugo anti-PD-1 (GSE78220), and Auslander (GSE115821). Cohort
selection for the ICB endpoint was gated on the availability of pre-treatment samples with
RECIST response labels; only Riaz was adequately powered (see Statistics).

## Single-cell analysis

### Single-cell processing and clone inference

Count matrices were loaded into AnnData and normalized to 10,000 counts per cell followed by
log1p (`code/load_adata.py`, `code/luad/build_luad.py`; for the 208,506-cell LUAD text
matrix a streaming line-by-line sparse reader was used to keep only tumor-derived cells with
adequate epithelial content, avoiding a dense in-memory matrix). Gene chromosomal positions
were taken from a BioMart GRCh37 export and attached to the variable index.

Malignant clones were inferred per patient with infercnvpy 0.6.1. For each patient, the
malignant compartment (epithelial or lineage-appropriate) was scored against a
patient-matched normal reference (non-malignant immune and stromal cells by author cell-type
annotation) using a 100-gene sliding window (`infercnvpy.tl.infercnv`). To run inside a
sandboxed environment the parallel map was monkeypatched to serial and the numba cache
directory was set before import. Cells were clustered on the CNV profile
(`infercnvpy.tl.pca` with `use_rep="cnv_pca"`, `pp.neighbors`, Leiden at resolution 0.4).
CNV clones were called malignant when they were majority-malignant by cell type, carried a
mean absolute CNV score above 1.5x the reference, and contained at least 30 cells. A
whole-tumor CNV ratio (mean malignant / mean reference CNV score) was recorded per patient.

### Non-circular validation

Because clones are defined from expression-derived CNV, we validated that the clone
partition reflects genuine transcriptional structure rather than a CNV-clustering artifact.
For each multi-clone tumor we re-embedded its malignant cells in an independent feature space
(2,000 highly variable genes, scaled, 20 principal components) and computed the silhouette of
the CNV-defined clone partition in that space. This was compared to a within-patient
permutation null (30 label permutations); we report the z-score of the observed silhouette
against the null (`code/phase0_validate.py`, `code/luad/run_luad.py`,
`code/aml/validate.py`).

### Clone pseudobulk, programs, and immune-state label

Per clone we formed a pseudobulk profile as the mean log-normalized expression across its
cells. Hallmark transcriptional programs (tissue-appropriate differentiation, proliferation,
EMT/invasion, an interferon/MHC antigen-presentation program, and cohort-specific oncogenic
programs such as the DDLPS MDM2/CDK4 12q amplicon) were scored per clone with
`scanpy.tl.score_genes` (`code/phase0_programs.py`, `code/luad/run_luad.py`). Recurrent
clone classes were obtained by Ward clustering of z-scored program scores at k=2. The
harmonized immune-state label ("immune-hot") was defined by orienting each cohort's classes
on the interferon/MHC program, so the label has consistent biological meaning across cohorts
whose raw class indices differ. For LUAD, whose interferon/MHC axis is continuous rather than
bimodal, the label was assigned by a median split.

### Effect-size and dilution analyses

To quantify imprint strength we compared within-patient subclone separation to
between-patient separation (mean absolute difference per gene across the top-2,000-variance
genes of the standardized pseudobulk matrix; `code/phase0_effectsize.py`). To assess
robustness to stromal/immune admixture we diluted malignant pseudobulk with reference
expression across a purity gradient and measured clone separability (silhouette) at each
purity (`code/phase0_dilution.py`). To test whether admixture can be countered we identified
malignant-compartment-specific genes (specificity = malignant vs reference expression
contrast) and repeated the low-purity separability test restricted to those features
(`code/phase0_compartment.py`, `code/ddlps_compartment.py`).

## Classification and transfer

### Supervised classifier

Clone-level pseudobulk profiles were pooled across cohorts on the shared gene set
(`code/pooled_classifier.py`). Features were z-scored within cohort to remove per-cohort
scale and batch offsets, then reduced to the top-2,000-variance genes. The classifier was L2
logistic regression (C=0.05, class-weight balanced) predicting the harmonized immune-hot
label. Generalization was estimated two ways: (i) grouped-by-patient stratified 5-fold
cross-validation, so no patient's clones appear in both train and test; and (ii)
leave-one-cohort-out, training on all tissues but one and testing on the held-out tissue, to
probe cross-tissue (zero-shot) transfer. Significance was assessed with a permutation null
(200 permutations) in which labels were shuffled within cohort and the full grouped-CV was
re-run each time.

### Real-bulk transfer and the ICB endpoint

On real DDLPS bulk (GSE221492) we scored the same programs and correlated them with the
matched-patient single-cell class, and separately tested NNLS-based composition
deconvolution of the bulk against a single-cell signature matrix
(`code/ddlps_pseudobulk.py`, `code/ddlps_compartment.py`). For the ICB endpoint
(`code/icb/riaz_icb_test.py`), Riaz FPKM (Entrez-indexed) was log2-transformed and
gene-z-scored across samples. A malignant-intrinsic interferon/MHC signature (HLA-A/B/C,
B2M, STAT1, IRF1, TAP1/2, GBP1, CXCL10, PSMB9, HLA-E, NLRC5) and an immune-infiltration proxy
(PTPRC, CD3D/E, CD8A, CD68, LYZ, MS4A1, NKG7, CD2, ITGAM, GZMB/A) were scored as mean z per
sample in the 49 pre-treatment samples. Each signature was tested against RECIST response
(responder = PRCR; non-responder = PD or SD) by Mann-Whitney U and ROC AUC, the two
signatures were correlated (Pearson), and a joint logistic model (response ~ infiltration +
IFN/MHC, statsmodels) tested whether the malignant-intrinsic signature adds beyond
infiltration.

### Statistics and power

A power analysis calibrated on real melanoma bulk covariance (Hugo GSE78220) established the
usable sample-size floor: a moderate within-sample effect needs roughly 40-67 matched
patients for balanced accuracy 0.62-0.82, and a strong effect is usable from about 25
(`code/` power grid). Reported classifier metrics are balanced accuracy and AUC under
grouped cross-validation with permutation-based p-values. Validation z-scores are relative to
within-patient label-permutation nulls. All random seeds were fixed (0) for reproducibility.

### Software

Python 3.11 with scanpy 1.11.5, infercnvpy 0.6.1, scikit-learn, statsmodels, scipy, and
anndata. All analysis scripts are under `analysis/clonal_trajectory/code/`, organized by
cohort (`aml/`, `luad/`, `icb/`) and by phase (`phase0_*`, `run_*`, `pooled_classifier.py`).
