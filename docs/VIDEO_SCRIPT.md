# Demo video — narration script + storyboard (≤ 3:00)

Target: **2:50**. Read at a calm pace (~140 wpm). Each shot lists what's on
screen and the exact words. Figures referenced are in `results/` and
`results/figure_deck.pdf`; the live demo is `notebooks/demo_reproduce_headline.ipynb`.

---

### Shot 1 — Title / hook (0:00–0:20, ~45 words)
**On screen:** `fig_hook_replication.png` (Fig 1).
> "Tumor-mutational burden and the other approved immunotherapy biomarkers top
> out around an AUROC of 0.6. The field keeps hoping RNA holds a better signal.
> We asked a sharper question: if you build an RNA biomarker of checkpoint
> response, does it actually *generalize* — or does it just look good on the
> cohort you trained it on?"

### Shot 2 — The thesis (0:20–0:40, ~45 words)
**On screen:** README thesis paragraph, or a simple slide of the hypothesis.
> "The hypothesis: early driver mutations set a tumor's evolutionary
> trajectory, and downstream RNA phenotypes — splicing, intron retention, RNA
> editing, transposable-element activation — are manifestations of one latent
> RNA state that shapes antigenicity and response. We tested it two ways, and
> tried hard to break each one."

### Shot 3 — Internal null (0:40–1:05, ~55 words)
**On screen:** `fig_covariation.png` (Fig 2) then `fig_power_curve` (Fig 3).
> "First, using standard exome-derived neoantigen proxies. If a shared RNA
> state existed, these should co-vary beyond mutational burden. They don't —
> permutation p of 0.78. And this wasn't underpowered: the test could detect a
> shared-variance of 23 percent; we observed 2. The proxies are the wrong
> instrument — exome annotations can't see a transcriptomic state."

### Shot 4 — Raw-read pivot (1:05–1:35, ~60 words)
**On screen:** `fig_incremental_auroc.png` (Fig 5), then a terminal snippet of
`run_salmon_pilot.sh` streaming.
> "So we went to the raw reads. No approved proxy adds anything over the TMB
> floor. We built a bandwidth-aware pipeline from scratch — salmon, arm64,
> streaming the first three million read-pairs per sample straight from ENA,
> quantifying, verifying, and deleting — and de-novo quantified 52 melanoma
> transcriptomes across three cohorts on a laptop."

### Shot 5 — The seductive signal (1:35–1:55, ~40 words)
**On screen:** `fig_pilot_denovo.png` (Fig 6).
> "And we found something. A de-novo antigen-presentation program — B2M, HLA,
> TAP, PSMB9 — cleanly separates responders within the Gide cohort:
> leave-one-out AUROC 0.87, robust to cross-validation scheme and to which
> genes we pick. Exactly the kind of result you'd want to believe."

### Shot 6 — The teardown (1:55–2:30, ~65 words)
**On screen:** `fig_infiltration_confound.png` (Fig 7), then
`fig_transfer_3cohort.png` (Fig 8).
> "Then we stress-tested it. First: that axis is 0.77 correlated with immune
> infiltration — every antigen gene tracks leukocyte content, so it's mostly
> reporting how much immune infiltrate is present. Second, the decisive test:
> trained on Gide, it scores 0.36 on held-out Riaz and 0.58 on held-out Hugo.
> A signal that looked real in one cohort vanishes in two others."

### Shot 7 — Payoff + reproducibility (2:30–2:50, ~40 words)
**On screen:** the 4-cell `demo_reproduce_headline.ipynb` running.
> "That's the honest finding: gene-level de-novo expression reflects the
> microenvironment and is cohort-specific — not a transferable RNA state.
> Testing the real hypothesis needs the WES-blind, multi-phenotype features
> standard pipelines discard, which we built and staged. Everything reproduces
> in a thirty-second notebook. Thanks for watching."

---

**Recording tips**
- Screen-record the figure deck PDF full-screen; advance one figure per shot.
- For Shot 4, a few seconds of the pipeline log scrolling sells the "from raw
  reads" claim.
- For Shot 7, run the notebook live — the four printed numbers land the arc.
- Total word count ≈ 390 words → ~2:45 at 140 wpm, leaving headroom under 3:00.
