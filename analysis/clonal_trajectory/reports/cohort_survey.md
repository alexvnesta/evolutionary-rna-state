# Cohort survey: closing the N≈40 matched-cohort gap

## The requirement

The investigation's single binding constraint is a **within-tissue cohort of ≥40 patients
with both scRNA-seq and matched bulk RNA-seq**. All three cohorts analyzed to date fell
short on the matched axis: cSCC 9 patients (WES cohort-pooled, no per-patient bulk),
DDLPS 19 matched, AML 16 (no bulk). The pooled N=99 result cleared the power floor only by
combining three tumor types, and leave-one-cohort-out showed that pooling does not buy
cross-tissue transfer. So the open question is specifically single-tissue power.

## What the GEO survey found

A structured search across melanoma, NSCLC, breast, ccRCC, colorectal, HNSCC, glioma, and
generic matched-atlas terms returned few series that carry both assays, and the strongest
single-tissue scRNA cohorts do not deposit matched bulk in the same series:

| Cohort | Tissue | scale | matched real bulk? |
|---|---|---|---|
| GSE131907 (Kim) | lung adenocarcinoma | 44 patients, 208k cells, 58 samples | not in series |
| GSE174554 | IDH-wt glioma | 80 snRNA specimens (40 primary + 40 matched recurrent) | not in series |
| GSE176078 (Wu) | breast | 26 patients | partial (CITE-seq, not bulk) |
| GSE115978 (Jerby-Arnon) | melanoma | 7186 cells | no (already used as ICB training ref) |

## Structural conclusion

Matched scRNA + real bulk from the same ≥40 patients is essentially absent from public
GEO as a single linkable deposit. Two realistic routes remain:
1. **Within-tissue powered same-cell test.** Use a large single-tissue scRNA cohort
   (GSE131907, 44 LUAD patients) and run the full clone-trajectory pipeline with
   within-cohort pseudobulk. This directly tests whether the immune-state axis is
   recoverable at N≥40 in ONE tissue, the claim currently only inferred from the pooled
   cross-tissue result. Reachable now with data in hand.
2. **Federated real-bulk assembly.** Cross-reference a scRNA cohort's patient IDs against
   TCGA / ICGC bulk for the same tumor type. This gives real bulk but breaks the
   same-patient link, so it can only test unsupervised structure and outcome correlation,
   not supervised label recovery.

Route 1 is the immediate next experiment; route 2 is the eventual deployment test.
