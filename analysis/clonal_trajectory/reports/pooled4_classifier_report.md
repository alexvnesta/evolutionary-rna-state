# Four-cohort pooled classifier: LUAD added

## What changed

The pooled powered classifier is re-run with lung adenocarcinoma (LUAD, GSE131907) added
as a fourth cohort, taking the pool from 99 to 145 clone-level pseudobulk
profiles (cSCC 25, DDLPS 55,
AML 19, LUAD 46). Design is
unchanged: harmonized immune-hot label oriented per-cohort by interferon/MHC, features
z-scored within cohort, top-2000-variance, L2 logistic, grouped-by-patient 5-fold CV plus
leave-one-cohort-out, permutation over labels within cohort.

## Result: more N, but a lower ceiling — tissue structure, not sample size, now binds

Grouped-by-patient CV gives balanced accuracy 0.79, AUC
0.87 (N=145, permutation p=0.005). This is
still strongly significant but lower than the three-cohort result (balanced accuracy
0.91, AUC 0.94 at N=99). Adding a fourth
tissue with more samples did not raise the pooled score; it lowered it.

The reason is documented in the LUAD replication: LUAD's interferon/MHC signal is
continuous rather than a discrete hot/cold split, so forcing it into a median-split label
injects noise the other three cohorts do not have. This is the key lesson: past the power
floor, the binding constraint stops being N and becomes how discretely each tissue's
immune axis is structured.

Leave-one-tissue-out (N=145): cSCC 0.57, DDLPS
0.61, AML 0.74, LUAD 0.61.
AML transfer improved markedly versus the three-cohort run (0.53 to
0.74) now that three tissues train it, but cSCC, DDLPS, and LUAD stay in
the 0.57-0.61 band. So more training tissues help a genomically-quiet target like AML but
do not deliver clean zero-shot transfer for the others: the method remains within-tissue.

## What this settles

- The three-cohort balanced-accuracy 0.91 is not a universal ceiling; it reflected three
  tumors that happen to share a discretely-structured immune axis. A fourth tumor with a
  graded axis pulls the pooled score down to 0.79 while staying highly significant.
- The practical design conclusion sharpens: train and apply within a tumor type, and
  expect the achievable accuracy to depend on how bimodal that tumor's immune-state
  structure is, not only on how many matched samples are available.

## Artifacts

- `figures/pooled4_classifier.png` : 3-vs-4-cohort CV, leave-one-tissue-out, permutation null.
- `tables/pooled4_results.json`, `pooled4_labels_preds.csv`.
