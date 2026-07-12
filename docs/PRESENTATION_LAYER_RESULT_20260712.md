# Presentation-layer campaign — presented novel-junction neopeptide load vs the immune floor

**Date:** 2026-07-12 · **Cohort:** Gide 2019 (within-cohort, n=30, 16R/14N) · **Status:** complete, null

## Motivation
Every non-reference RNA feature we tested (splicing, editing, IR, TE, Evo 2 sequence
aberrancy, aberrancy×dosage interaction) was null beyond the immune floor, and the unifying
mechanism was that bulk "aberrancy" is entangled with immune composition and inherits the
cross-cohort sign-flip that sinks the floor. The one untested salvage was the **presentation
layer**: aberrant neopeptides filtered through (a) NMD-escape and (b) the patient's own HLA —
a feature that should, in principle, be **decoupled from bulk infiltration**.

## Pipeline (built this session, parameterized by cohort)
1. **Novel-junction neopeptides.** For each of the 1,864 top-200 novel junctions per Gide
   sample, translate the spliced sequence in all 3 frames and collect 8–11-mer peptides that
   **span the splice site**.
2. **NMD-escape filter.** Retain only peptides from ORFs that plausibly escape nonsense-
   mediated decay (premature stop codon read-through, or stop within the EJC-proximal 55 nt of
   the modeled acceptor window). This removes ~92% of naive junction peptides
   (20,269 → 1,714 unique 9-mers; 429 of 1,864 junctions retain ≥1 escaping peptide).
   *Limitation:* first-order model using the 512 bp acceptor flank as an EJC proxy, not a full
   transcript reconstruction.
3. **Patient HLA.** arcasHLA typing of all 30 local Gide markdup BAMs (HLA-A/B/C, 4-digit).
   Mean class-I heterozygosity fraction 0.822. Table: `results/hla_typing/gide30_hla_alleles.csv`.
4. **MHC-I presentation.** mhcflurry `Class1PresentationPredictor` scores each patient's
   NMD-escaping neopeptides against their own alleles → per-sample presented-load features
   (count of strong presenters, read-weighted presented load, etc.).
5. **Two-block test.** Same rigorous protocol as the encoder tests: 20-seed 5-fold CV,
   fold-contained residualization, cohort-internal permutation null; plus a rich-basis
   disentangling frame (floor + 11 InstaPrism deconvolution fractions).

## Result

| Metric | AUROC |
|---|---|
| Immune floor (positive control) | **0.854** |
| Presented load alone | 0.581 |
| Floor + presented | 0.846 (Δ −0.008) |
| Presented residualized on floor (fold-contained) | 0.472 |
| — cohort-internal permutation p | **0.552 (NS)** |
| Presented residualized on rich immune basis (floor + deconv) | 0.413 |

**Verdict: null.** Presented-neopeptide load carries no independent ICB-response signal beyond
the immune floor within Gide. Residualized it sits at/below chance (0.472, p=0.55) and drops
further against the rich deconvolution basis (0.413) — the same proxy signature as every other
non-reference feature. The strongest formulation of the hypothesis does not rescue it.

## Scope and limits
- **Single cohort, n=30.** Hugo has zero local BAMs, so no within-cohort replication and no
  transfer test. This is a within-Gide result, consistent with (not independent of) the prior
  two-cohort encoder nulls.
- A null at n=30 bounds but does not exclude a small effect a larger cohort might reveal — which
  is exactly what the in-progress 106-sample alignment on `ubuntuserver.local` will enable. The
  pipeline is parameterized by cohort and fires unchanged at 106-scale.
- The NMD model is first-order (see step 2).

## Artifacts
- `analysis/build_presented_load.py`, `analysis/presented_two_block.py`, `pipelines/scripts/hla_type_gide.sh`
- `results/predictor/presented_block_gide2019.parquet`
- `results/eval/presented_two_block_gide2019.json`, `results/eval/presented_perm_null_gide2019.npy`
- `results/hla_typing/gide30_hla_alleles.csv` / `.parquet`
- `presented_layer_result.png`
