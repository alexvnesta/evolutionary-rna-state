# Results — viral-mimicry structure & latent-evolutionary-state probes

**Date:** 2026-07-22 · **Session:** `c55b6dde` · **Motivated by:** `GEMINI_ASSESSMENT_20260722.md`
**Verdict in one line:** Both prongs are **negative** — but each negative is *sharper and more
informative* than the burden-level nulls that preceded it. The dsRNA structure feature is the
first non-reference feature that is genuinely **decoupled from the immune floor** yet still
carries **no ICB signal**; and the RNA latent **does not reconstruct the tumor's clonal
architecture** in the first properly-powered test of the hypothesis's literal claim.

Figure: `results/fig_gemini_twoprong.png`

---

## Prong A — inverted-repeat dsRNA / viral-mimicry ΔG feature

**Rationale.** We measure A-to-I editing and TE/ERV family *abundance*, but never the physical
quantity that triggers innate dsRNA sensing (RIG-I/MDA5/TLR3): whether expressed retroelement
transcripts fold into stable double-stranded structure. Gemini named this ("inverted-repeat
pairing energy"); it is the one genuinely new, mechanistically-motivated feature in that exchange.

**Construction.**
- Genome-wide, found proximal **opposite-strand Alu pairs** (gap ≤ 1000 nt, both arms inside one
  expressed gene body) — opposite-strand copies of the same repeat are reverse-complementary and
  form dsRNA hairpins (the IRAlu model).
- Scored each pair's hybridization free energy with **ViennaRNA `duplexfold` (ΔG)**:
  **286,368 pairs across 20,380 genes.** Per-gene: Σ/mean/min ΔG, pair count.
- Per-sample **dsRNA structure signature** = mean z-scored log-expression over the top-decile
  IR-Alu-structure gene set (2,038 genes), then **residualized on global expression in-fold** so
  the feature is decoupled from overall expression magnitude (corr with total log-expr → 0.0).
  This is a *structure* feature, not a repeat-abundance count: the genomic IR-Alu map is fixed
  across samples; only expression weighting varies.

**Decoupling diagnostic (the decisive question).** A viral-mimicry index is *expected* to track
interferon (it signals through the IFN axis = our immune floor). It does **not**:

| dsRNA structure vs | Gide (30) | Hugo (25) | Riaz (10) | All (68) |
|---|---|---|---|---|
| IFN-γ score (ρ) | −0.045 | −0.248 | −0.139 | −0.090 |
| GEP T-cell-inflamed (ρ) | −0.088 | −0.178 | +0.042 | −0.080 |

All |ρ| < 0.26, none significant. **This is the first non-reference feature in the project that is
not an inflammation proxy** — every prior non-ref layer (aberrancy burden, novel junctions,
editing) correlated with and sign-flipped against the immune axis (`EVO2_INTERACTION_INSIGHT`).

**Two-block test (20-seed 5-fold CV, fold-contained residualization):**

| | Gide n=30 | Hugo n=25 |
|---|---|---|
| immune floor | 0.853 ± 0.014 | 0.677 ± 0.039 |
| dsRNA structure alone | 0.310 ± 0.050 | 0.331 ± 0.074 |
| floor + dsRNA | 0.854 (Δ +0.001) | 0.638 (Δ −0.039) |

(Both blocks residualize the dsRNA feature on global expression **fold-contained** — train-fit only.)
Residualized-on-floor + cohort-internal permutation (500×): Gide AUROC 0.400 (perm p=0.643),
Hugo 0.514 (perm p=0.363) — both at/below the null. The below-chance "alone" values are small-n
noise, confirmed by permutation.

**Interpretation.** The dsRNA structural feature is **decoupled-but-null**: it is not the immune
floor in disguise, yet it carries no independent response information and adds nothing to the
floor. This is a *stronger* negative than a burden count — it says the strongest biophysical form
of the viral-mimicry hypothesis, built as a genuine structure feature orthogonal to inflammation,
still does not move ICB prediction in bulk RNA at this n. (`results/eval/dsrna_viral_mimicry_test.json`)

---

## Prong B — latent evolutionary state vs clonal ground truth

**Rationale.** The project hypothesis's literal claim is that bulk RNA contains enough information
to *reconstruct the tumor's latent evolutionary state.* The concrete evolutionary ground truth is
**clonal architecture** (from WES). `FINAL_VERDICT.md` flagged this was only ever probed at n=10
(all |ρ|<0.19, underpowered). **This session found it is testable at n≈60**, routing through
expression/scVI coverage rather than local raw-read quant only.

**Design.** Purified toward the tumor-intrinsic signal (InstaPrism **Malignant** fraction as a
covariate, n=106), RNA representation = scVI 30-dim latent (primary) + HVG-2000 expression PCA-20
(robustness). Targets = `heterogeneity_index`, `subclonal_fraction`, `mean_ccf`.
**Within-cohort only** — Hugo (purity-corrected CCF) and Riaz (purity-free VAF proxy) clonality are
computed by non-comparable methods (`_joins_provenance.json`); Gide has no WES.

**Result — within-cohort CV (multi-seed 5-fold RidgeCV, OOF Spearman, 500× permutation):**

| target | Hugo n=27 (obs / null / p) | Riaz n=33 (obs / null / p) |
|---|---|---|
| heterogeneity_index | −0.270 / −0.225 / 0.489 | −0.230 / −0.182 / 0.535 |
| subclonal_fraction | −0.316 / −0.210 / 0.631 | −0.422 / −0.196 / 0.932 |
| mean_ccf | −0.298 / −0.202 / 0.611 | −0.391 / −0.193 / 0.876 |

Every observed OOF ρ **equals its permutation null** (all perm p 0.49–0.93). The strongly-negative
nulls are the known small-n regularized-CV artifact; the point is observed ≈ null everywhere.
Best-single-dim in-sample upper bound (optimistic, ignores selection): Hugo tops at |ρ|≈0.41 (at
the multiple-comparison null for n=27); the higher Riaz expression-PC values are on n=10 coverage
only. Nothing survives.

**Interpretation.** **NULL, and now properly powered.** The RNA representation does not reconstruct
the tumor's clonal/evolutionary state at achievable within-cohort n (21–33). This is the first
powered test of the hypothesis's literal "reconstruct the evolutionary state" claim — and it is
tumor-intrinsic by construction (purified, anchored to WES, unsupervised representation, no response
label doing the fitting), so it is *not* confounded by immune composition the way every
response-anchored test was. It simply is not there in bulk RNA at this scale.
(`results/eval/latent_state_clonality_probe.json`)

---

## What this settles, and what it does not

- **Settles:** Gemini's one novel idea (dsRNA structure) and Alex's reframe (model the latent state,
  anchored to clonality) are the two most promising untested directions — and both, tested cleanly,
  are negative. The viral-mimicry feature is not a hidden immune proxy; it is genuinely absent as an
  independent signal. The latent-state claim does not hold at n≈60 tumor-intrinsic.
- **Does not settle:** (i) A cross-cohort clonality test at full power needs Hugo+Riaz recomputed
  with one pipeline (PyClone-VI + FACETS/ASCAT purity) and ideally WES on Gide (absent). (ii) The
  dsRNA feature is Alu-IR only; ERV/LINE dsRNA and editing-modulated stability (ADAR destabilizes
  dsRNA) are not yet in it. (iii) Single-cell would resolve the malignant compartment far better
  than bulk deconvolution. None of these is likely to overturn the direction given the consistency
  of the nulls, but they bound the claim.

**Bottom line.** Two sharp, falsification-first negatives. Consistent with the project's converged
position: the ICB signal in bulk tumor RNA is immune composition; tumor-intrinsic non-reference and
evolutionary-state features — even in their strongest, best-motivated, confound-decoupled forms —
do not carry independent response information at achievable n.
