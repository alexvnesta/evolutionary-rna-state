# LUAD (Kim GSE131907) — fourth-tissue replication and largest single cohort

## Why this cohort

The investigation's one binding constraint is a large single-tissue cohort. A GEO survey
(cohort_survey.md) found that matched scRNA + real bulk from the same >=40 patients is
essentially absent from public repositories, but it identified GSE131907 (Kim lung
adenocarcinoma, 208,506 cells) as the largest single-tissue scRNA cohort with author
cell-type annotation. Restricting to tumor-derived samples with at least 200 epithelial
cells gives 26 tumors, the largest single cohort in this project (cSCC 9, DDLPS 19,
AML 16). Lung adenocarcinoma is also a fourth, glandular-epithelial tumor type distinct
from squamous cSCC, mesenchymal DDLPS, and hematopoietic AML.

## Result 1: the clone-imprint replicates, strongly

Of 26 tumors, 25 show a whole-tumor CNV ratio at or above 1.5x (median
2.25, range 1.24-3.42);
LUAD is genomically active, unlike AML. Clone inference found 16 multi-clone
tumors. The non-circular validation (re-clustering malignant cells in independent
highly-variable-gene expression space, scored against a within-patient permutation null)
puts 15 of 16 multi-clone tumors above the null (median z 9.9, up to
134). Within-patient subclone separation is 0.78 SD/gene versus
1.05 between patients, the same pattern as the other three cohorts.

## Result 2: the immune axis is present but continuous, not a discrete hot/cold split

This is where LUAD differs. In cSCC, DDLPS, and AML the clones organized into two discrete
recurrent classes with an immune-hot pole carrying roughly 36% of clones. In LUAD the
interferon/MHC signal across the 46 clones is unimodal and continuous (no clean gap; a
blind two-class clustering isolates a small alveolar-differentiated outlier rather than an
immune dichotomy). The oncogenic programs do not anti-correlate with immune signaling the
way the DDLPS MDM2 amplicon did (r=-0.63 there); here proliferation shows only a weak
negative trend (r=-0.28, p=0.06) and EMT is weakly positive (r=+0.38, p=0.01). So the
two-state immune axis is tumor-type-dependent: sharp in cSCC/DDLPS/AML, graded in LUAD.

## What this adds

- A fourth independent replication of the core clone-imprint result, in the largest single
  cohort assembled here and in a fourth tumor lineage. The imprint itself (clones are real,
  genome-wide, non-circular, recurrent) is now robust across four cancers.
- A boundary condition on the immune-class structure: the discrete immune-hot/cold
  dichotomy is not universal. LUAD's antigen-presentation variation is continuous. This
  matters for the pooled classifier, whose harmonized label assumes a hot/cold split; LUAD
  is included via a median split, and its softer axis is expected to lower cross-tissue
  transfer to and from LUAD specifically.

## Artifacts

- `figures/luad_replication.png` : CNV ratio, independent validation z, IFN/MHC distribution.
- `tables/luad_infercnv_summary.csv`, `luad_validation.csv`, `luad_clone_programs_classified.csv`, `luad_effectsize.csv`.
- `tables/luad_clone_pseudobulk.npy` (+meta, genes) for the pooled classifier.
- `code/luad/build_luad.py`, `code/luad/run_luad.py`.
