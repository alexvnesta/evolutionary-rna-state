# Assessment of the Gemini advice exchange — and the plan it motivates

**Date:** 2026-07-22 · **Author session:** `c55b6dde` · **Status:** assessment + execution plan
**Trigger:** Alex shared a multi-turn Gemini conversation (pasted-text 2026-07-22) that
began with a "what is MLP" question against a live agent log and evolved into a design
discussion about combining self-supervised learning, pretrained sequence LMs, and feature
engineering to model the tumor "latent state." This document records what in that exchange
is correct, what is misdiagnosed, and the two concrete experiments it motivates — both of
which are runnable on present local substrate.

---

## 1. Citation audit (done first, because specific paper+year claims are where chatbots fabricate)

Every load-bearing citation Gemini gave was checked against a retrieved source. They hold up:

| Claim | Status | Source |
|---|---|---|
| P-NET — biologically informed sparse NN, prostate cancer | **CONFIRMED** | Elmarakeby et al., *Nature* 598:348–352 (2021), doi:10.1038/s41586-021-03922-4 |
| LEMBAS / BINN robust to spurious interactions via self-pruning | **CONFIRMED** | bioRxiv 2025.10.24.684155 (Oct 2025) |
| Pathway-informed VAE increases interpretability | **CONFIRMED** | Liu et al., *PLOS Comput Biol* 20:e1011198 (2024) |
| bMAE — masked autoencoder for bulk RNA-seq tissues (GTEx pretrain) | **CONFIRMED** | bioRxiv 2026.03.03.709470 |
| "TxFM — masking gene expression" | **UNCONFIRMED** | no matching record surfaced; treat as unverified |
| Viral mimicry / HERV / dsRNA sensing → ICB | **CONFIRMED** (established biology) | Cancer Discovery 11:2707; JCI 183745; MDPI Biomolecules 16:709 |

So the literature scaffolding is sound. The problem is not fabricated sources — it is that the
advice is aimed at the wrong failure mode. That is the substance of this assessment.

## 2. Where Gemini is right (and consistent with our own results)

- **"The aberration block hurts, linear or nonlinear."** Gemini's own MLP run reproduced our
  central null with a fresh model class — pooled AUROC 0.61→0.46 when the 35 aberration axes are
  added. This matches the documented non-ref-degrades-floor result (within-Gide Δ −0.271;
  the Evo2 / EVA / HyenaDNA encoder nulls; `DISENTANGLE_TEST_20260711.md`). Corroboration, not novelty.
- **"Do not fine-tune a 7B sequence LM end-to-end on n≈100."** Correct, and already our position
  (`ENCODER_PHASE_PROTOCOL.md`: frozen embeddings first, PEFT gated, training signal aimed at
  phenotype reconstruction, never directly at response).
- **"Modular / multi-tier, not single end-to-end."** Correct, and again the protocol's design.
- **The gradient-attenuation / scale-mismatch / leakage argument** against a single end-to-end
  net is a fair, if florid, statement of the small-n overfitting trap.

None of this is wrong. Most of it we already knew and wrote down.

## 3. Where the diagnosis is importantly wrong — the wall is confounding, not overfitting

Gemini's framing, repeated throughout, is that the wall is **dimensionality / variance at n≈92**
("35 axes overwhelm the sample size," "curse of dimensionality," "the model learns nonlinear
noise"). Almost all its prescriptions follow from that premise: masked networks as regularizers,
aggregation to 2–3 scalars, SSL pretraining to smooth variance.

Our own evidence says that is **not** the primary failure mode. Three results Gemini did not have
point elsewhere:

1. **Aberrancy×dosage sign-flip** (`EVO2_INTERACTION_INSIGHT_20260711.md`): the non-ref burden
   correlates with every inflammation axis, and its coupling to inflammation *flips sign across
   cohorts* (GEP +0.36 in Gide → −0.15 in Hugo). This is a transportability failure, not a
   variance one.
2. **Disentangle test** (`DISENTANGLE_TEST_20260711.md`): residualizing the aberration signal on
   the rich 16-dim InstaPrism immune basis, fold-contained, leaves an orthogonal subspace with
   **no response information** — it drops toward chance rather than being rescued.
3. **Purity-corrected mediation** (`AUDIT_SYNTHESIS.md` §Phase 1b): the phenotype→infiltration
   a-path collapses to ~0 once tumor purity is removed.

Together these say the aberration feature is a genuine **downstream proxy of the tumor-immune
state whose relationship to that state reorients between cohorts.** That is a confounding +
non-transportability problem. It has two consequences for Gemini's advice:

- **Its two top fixes are experiments we already ran, negative.** "Regress out immune signatures,
  feed the residual" and "orthogonal projection" *are* the disentangle test. The residual subspace
  is empty. A regularizer cannot manufacture signal where the orthogonal subspace carries none.
- **A biologically-masked VAE does not escape the confound — it re-learns it.** Routing aberration
  features into a "viral mimicry node" will simply have that node covary with the immune node,
  because that is the true biology (viral mimicry acts *through* interferon). Masking fixes
  spurious cross-talk; it does not create an independent axis where the data has none.

And **more scale does not touch transportability.** The CheckMate-064 push to n≈128 buys power —
worth having — but four entangled cohorts test the sign-flip more stringently; they do not remove
it. Do not let the n≈128 alignment reframe the question as a power problem.

## 4. The one genuinely valuable new idea — structure-aware dsRNA / viral-mimicry ΔG

The **inverted-repeat complementarity / dsRNA fold-stability (ΔG)** feature is the single
suggestion that is both new relative to our pipeline and mechanistically motivated. We measure
A-to-I editing (AEI) and TE/ERV *family abundance*, but we have never computed the physical
quantity that actually triggers innate sensing: whether expressed retroelement transcripts form
stable double-stranded structure (RIG-I / MDA5 / TLR3 ligands). That is a different feature *in
kind* from a bulk burden.

**The caveat Gemini did not flag, and it is decisive:** viral mimicry signals *through* the
IFN / dsRNA-sensing axis → antigen presentation and infiltration — and that IFN / T-cell-inflamed
axis **is our immune floor.** So a dsRNA index is mechanistically expected to correlate with the
very floor every non-ref feature has failed to beat. It may well be one more confounded proxy — a
better-motivated one. It is worth building *because it is a sharp, cheap, falsifiable test of the
strongest biophysical form of the hypothesis*, not because it is likely to be the exception. A
null that collapses onto the floor is a **stronger** negative than a burden count.

→ **Prong A** below.

## 5. The most important move — Alex's reframe toward "model the latent state," done correctly

Alex's pivot — *"I don't want to just predict ICB response, I want to model the latent state of
these tumor samples"* — is the right instinct and sidesteps the trap. Everything to date regresses
features against **response**, which is dominated by immune composition, so every route
rediscovers composition. An *unsupervised* representation anchored to something tumor-intrinsic is
not fighting that confound.

**But the caveat neither Gemini nor the thread confronts:** an unsupervised manifold of *bulk*
tumor RNA is *also* composition-dominated — atlas-scale bulk VAEs organize primarily by tissue and
composition (bMAE, the pathway-VAE both show this). So "model the latent state" on bulk RNA, done
naively, gives back a composition manifold — the same trap in unsupervised clothing. The
evolutionary state we want is a small tumor-intrinsic residual *underneath* composition. The lever
that surfaces it is not architecture — it is **what you anchor to and what you purify to:**

1. **Anchor to a real evolutionary ground truth, not response.** The literal latent variable is
   clonal architecture. `FINAL_VERDICT.md` flags this was only ever probed at n=10 (RNA-axis vs
   clonality, all |ρ|<0.19, underpowered). **On-disk feasibility recheck (this session): it is not
   n=10 — it is n≈60.** Routing through expression/scVI coverage (not local raw-read quant only),
   the RNA∩clonality set is Hugo 27 + Riaz 33 = 60 samples with WES-derived clonality
   (`hugo_clonality.csv`, `riaz_clonality.csv`: heterogeneity_index, subclonal_fraction, mean_ccf).
2. **Purify before you represent.** Use the InstaPrism **Malignant** fraction (already computed,
   n=106, `instaprism_fractions_n106.parquet`) so composition is not the dominant axis.

A masked / pathway-informed VAE (Gemini's P-NET / pathway-VAE direction) is *reasonable here* in a
way it is not for response prediction: with reconstruction + evolutionary-anchor objectives and no
confounded response label doing the overfitting, the latent state is defined by tumor-intrinsic
structure. Build it to reconstruct/anchor the **evolutionary** state on the **malignant
compartment** — not to predict response.

**Hard honesty constraint on Prong B:** Hugo and Riaz clonality are computed by *non-comparable*
methods (Hugo purity-corrected CCF; Riaz purity-free VAF/clonal-peak proxy — see
`_joins_provenance.json`). Gide has no WES at all. So n≈60 is **two within-cohort probes**
(Hugo n≈21 usable after low-purity exclusion; Riaz n≈33), not one pooled test. Report effect sizes
+ CIs within cohort; no cross-cohort clonality pooling until both are recomputed with one pipeline.

→ **Prong B** below.

## 6. The plan this motivates

- **Prong A — viral-mimicry structure.** Build the inverted-repeat dsRNA/ΔG index from expressed
  TE/repeat loci (ViennaRNA on Alu/rmsk pairs, expression-weighted). Two-block test vs the immune
  floor (fold-contained, 20-seed CV, cohort-internal permutation) **plus** the decisive diagnostic:
  its correlation with the IFN/GEP axis. Expect entanglement; a null is a strong result.
- **Prong B — latent evolutionary state.** Purify to the malignant compartment (InstaPrism), build
  the RNA representation, and test whether it reconstructs the tumor's evolutionary ground truth
  (heterogeneity_index / subclonal_fraction / mean_ccf) *within cohort* at n≈60. This is the
  hypothesis's literal "reconstruct the evolutionary state" claim — tumor-intrinsic by construction
  and never powered before.
- **Synthesis.** Dated results doc + figures + PROJECT_STATUS milestone, honest verdict either way.

**Bottom line.** Gemini delivered a well-cited, competent tour of biologically-informed
architectures and independently reproduced our null — but it is solving overfitting, and our
problem is confounding. Alex's reframe toward modeling the latent state is the more promising path,
*provided it is anchored to a tumor-intrinsic ground truth rather than allowed to re-derive
composition.* The two prongs test exactly the two things neither the interpretable nor the encoder
work has yet addressed.
