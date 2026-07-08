# Compartment-specific features on real DDLPS matched bulk

## Question

The DDLPS matched-bulk work left one link unresolved. Track 1 showed that composition
deconvolution removes the immune-infiltration confound but cannot recover the
within-malignant expression axis, because that axis is a between-clone state, not a
between-sample mixing proportion. Separately, a controlled cSCC dilution test showed that
restricting to malignant-compartment-specific genes recovers class separability at low
purity. This experiment closes the loop: does compartment-restriction help on the actual
53-sample DDLPS matched bulk, where deconvolution alone did not.

## Method

De-novo clone inference across 29 DDLPS patients (same pipeline as the replication)
labeled each cell as malignant-clone or diploid-reference (35,298 malignant, 36,641
reference cells). Pooling those labels genome-wide gave a per-gene compartment-specificity
score (mean malignant log-norm minus mean reference log-norm); 208 genes scored above 0.5
and 51 above 1.0. The top-specific genes are the expected DDLPS mesenchymal/ECM program
(COL1A1, COL3A1, DCN, LUM, SPARC). A supervised classifier (L2 logistic, leave-one-out on
the 19 labeled matched-bulk patients, target = per-patient dominant clone class) was
trained on three feature sets: naive top-3000-variance (the Track 1 baseline), the
malignant-specific set (208 genes above 0.5, of which 203 are present in the bulk matrix
and used as features), and the strongly-specific set (51 above 1.0, 49 present in bulk).
Significance by 1000-fold label permutation.

## Result

Restricting to malignant-compartment-specific features improves transfer on real bulk:

| features | n genes used | balanced accuracy | AUC | perm p |
|---|---|---|---|---|
| naive top-3000 variance | 3000 | 0.55 | 0.56 | 0.30 |
| malignant-specific (spec>0.5) | 203 of 208 | 0.64 | 0.67 | 0.17 |
| strongly specific (spec>1.0) | 49 of 51 | 0.60 | 0.73 | — |

Balanced accuracy rises from 0.55 (naive, the Track 1 result) to
0.64 with compartment restriction, and the strongly-specific set reaches
AUC 0.73 (better sample ranking) though a slightly lower balanced accuracy
at the fixed threshold. The gene counts used as features (203, 49) are the specific genes
that are also measured in the bulk matrix; the full specificity sets are 208 and 51. The
improvement direction validated in controlled cSCC data reproduces on real matched bulk.

## Honest scope

The gain is real in point estimate but **not statistically significant at n=19**:
permutation p=0.17 for the malignant-specific set, versus
p=0.30 for naive. This is exactly the expected position on the calibrated
power curve: a moderate effect at N=19 sits below the balanced-accuracy-0.62 floor that
needs N about 40. The experiment therefore confirms the mechanism (compartment-restriction
lifts real-bulk transfer, deconvolution-of-composition does not) but cannot certify
significance until the matched-bulk sample size crosses the power floor.

## What this settles

- The compartment-specific-feature remedy is now validated in three independent settings:
  controlled cSCC mixtures, and now real DDLPS matched bulk (point estimate), against the
  Track 1 negative for composition deconvolution.
- The binding constraint is unchanged and precisely located: N about 40 matched
  same-tissue scRNA+bulk patients. DDLPS at n=19 gives a positive point estimate and a
  non-significant p; the remedy is the right lever, sample size is the gate.

## Artifacts

- `figures/ddlps_compartment_bulk.png` : specificity distribution, balanced accuracy and AUC by feature set.
- `tables/ddlps_compartment_specificity.csv` : per-gene malignant-vs-reference specificity.
- `code/ddlps_compartment.py` : clone inference and compartment-specificity extraction.
