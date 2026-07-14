# Non-reference ICB predictor — reasoning & open questions (session 15defe54, 2026-07-10)

_Durable notes for future agents. Captures the interpretation of the campaign's
results and the design decisions made in this session's Q&A. Numbers are quoted
from saved artifacts (loco_results.json, loco_results_scvi.json, imvigor_h2h.json)._

## 1. What "the predictor works at a modest immune-composition ceiling" means
- A predictor built purely from bulk RNA-seq DOES work, at a modest ceiling:
  immune-composition scores (Ayers T-cell-inflamed GEP + IFN-gamma + Mariathasan
  Teff/TGF-beta) reach AUROC ~0.66 on IMvigor210 (within-cohort CV) and ~0.54
  cross-cohort in melanoma LOCO (perm p=0.005 at n=136). This is a real,
  raw-RNA-derived predictor. The ~0.6-0.7 ceiling is BIOLOGICAL, not a coding
  limitation — even commercial models (GEM-1) live in that range on these cohorts.
- The tumor-INTRINSIC expression representation is at chance cross-cohort
  (linear latent 0.485 p=0.48; scVI 30-d 0.489 p=0.29). Same verdict from a
  proper neural VAE as from linear PCA — the null is robust.

## 2. Can we "correct for" immune composition? — THREE distinct operations
1. **To improve prediction: NO.** Immune composition is not a confound on top of
   the signal — as far as this data shows, it IS the signal. Regressing it out
   removes the predictive component and pushes AUROC toward chance. The ceiling
   is not contamination.
2. **To test incremental value of the tumor block: this is the RIGHT framing, and
   it is already built** as the `floor+X` vs `floor` increment (not subtraction —
   put the floor IN the model, ask if X raises AUROC on top).
3. **To rule out immune-composition as a hidden driver of the non-ref features:
   the composition-confound VALIDITY CHECK** — see section 4.

## 3. What "floor + nonref vs floor" (incremental test) has taught us so far
The increment is NULL-TO-NEGATIVE — adding a tumor block never significantly
improves on the immune floor, and tends to DILUTE it:
- immune_floor alone:            0.538  p=0.005   (cleanest signal in the campaign)
- floor + expression latent:     0.546 (scVI) / 0.561 (linear), p=0.095 / 0.134
  -> point estimate nudges up but the p-value gets WORSE — the latent adds noise
     faster than signal, eroding significance.
- floor + nonref:                NOT yet run at power (nonref was n=16; n=40
  preview gave 0.336 on a broken Hugo-heavy convenience subset — uninterpretable).
LESSON: the immune floor is a tight signal; everything tumor-intrinsic bolted on
so far has diluted rather than augmented it.

## 4. Composition-confound validity check (BUILT this session, twoblock_loco.py v2)
The sharp worry: the aberrant-RNA features may be DRIVEN BY immune infiltrate,
not tumor biology. ADAR is interferon-inducible, so an inflamed tumor shows more
A-to-I editing; infiltrate also shifts apparent TE/splicing. So if the non-ref
block ever shows signal, we must ask: is it tumor evolutionary state, or just
immune content re-detected through a different lens?
- IMPLEMENTED: `loco_auroc(..., resid_basis=floor)` residualizes each non-ref
  feature against the immune-floor scores PER-FOLD (linear fit on TRAIN only),
  then re-tests. Emitted as result key `nonref_resid_of_floor`.
- INTERPRETATION:
  - raw nonref_block significant BUT nonref_resid_of_floor -> chance
    => the signal WAS composition, NOT tumor-intrinsic (hypothesis NOT supported).
  - nonref_resid_of_floor stays significant
    => genuinely tumor-intrinsic, independent of composition (hypothesis supported).
- This runs AUTOMATICALLY alongside the incremental test when the powered matrix lands.

## 5. Labels/endpoints tested — and the GAP
- TESTED: binary RECIST response only (CR/PR = responder vs SD/PD). Two settings:
  melanoma anti-PD-1 3-cohort LOCO (Gide/Hugo/Riaz), bladder anti-PD-L1 (IMvigor210).
- NOT TESTED: any survival endpoint. **results/analysis_frame.parquet HAS
  OS_MONTHS/OS_STATUS/PFS_MONTHS/PFS_STATUS** — immune signatures often predict
  OS/PFS more sensitively than RECIST. Testing a Cox/survival endpoint against the
  floor and the non-ref block is a genuine untested angle, cheap to add.
- Also untested: durable clinical benefit (e.g. 6-month PFS), continuous BOR.

## 6. "Is there no way to test raw-RNA + foundation model?" — clarification
- Q1 "predict response from bulk RNA-seq at all?" -> YES, done, works at ~0.66 (floor).
- Q2 "does an existing FM embedding ADD anything?" -> tested, so far NO. CRITICAL:
  no FM ingests raw reads — Orthrus=mature-transcript SEQUENCE, Evo2=genomic DNA,
  scGPT/scVI=expression vector. You always quantify first. The scVI latent (the
  open GEM-1 analog) is at chance cross-cohort and does not beat the floor.
- Q3 the project's ACTUAL distinctive hypothesis (aberrant/non-reference RNA
  features carry signal) is GENUINELY UNTESTED at power — every negative so far
  rested on expression or a DNA antigen-quantity proxy, never the non-ref features.
  This is a LOGISTICS blocker (canonical per-feature cohort still processing), not
  a fundamental one.

## 7. The one untried route that could actually move the needle
Nobody has fed the ABERRANT TRANSCRIPT SEQUENCES THEMSELVES through a sequence
model. Everything tested treats the non-ref features as scalar summaries
(AEI %, TE CPM, IR fraction). The hypothesis's logic is that these sequences ARE
the antigen source — so embedding retained introns / chimeric-fusion transcripts /
edited-Alu elements as SEQUENCE (via Orthrus/HyenaDNA) and testing those
embeddings against response is conceptually distinct from everything done so far.
Speculative; needs the per-feature pipeline output first. This is the honest
answer to "is there anything left that could work."

## 8. Current honest status (calibrated)
- Prior on the hypothesis is NOT good: repeated proxy-falsifications; first weak
  look at the non-ref features (n=40, deprecated -k10, Hugo-heavy) shows nothing,
  though it is uninterpretable because the immune floor ALSO collapsed there
  (0.356) — the positive control failing means the subset is broken, not the biology.
- Status = "unpromising prior, but genuinely UNTESTED at power" — NOT "promising",
  NOT "falsified".
- What would change it: the FAIR test — canonical per-feature features (unique-mapper
  editing + bowtie2 -k100 TE-locus, where the hypothesis actually lives; family-level
  featureCounts discards the locus resolution that matters) on a Gide-representative
  cohort, with the immune floor recovered above ~0.55 as a validity gate, PLUS the
  composition-residualization check from section 4.

## Harness (all saved as artifacts, restore from store after workspace cleanup)
- assemble_nonref_block.py  (81bce881) — per-sample callers -> normalized block
- twoblock_loco.py          (5f6db3fd, v2) — leakage-guarded LOCO + perm + composition check
- nonref_attribution.py     (99ea392f) — per-feature mechanism attribution
- preview_loco_n40.json     (56c7faf2) — the n=40 preview (uninterpretable, documented why)
- BUILD_REPORT.md           (4141516d, v5) — full build report
RERUN TRIGGER: full 3-layer coverage crosses ~80 samples AND immune floor recovers
>~0.55 on the accumulated subset -> rerun twoblock_loco.py on the canonical block.


## 9. NEW FINDING (2026-07-10): fusion-neoantigen load associates with WORSE OS
Ran the survival endpoints that were never tested (section 5 gap). Data:
`results/analysis_frame.parquet` DNA/WES-derived non-reference neoantigen counts
(SPLICE/ERV/FUSION_NEOANTIGEN) + OS/PFS, n=264 across 4 ICB cohorts
(liu2019 122 / gide2019 72 / riaz2017 44 / hugo2016 26). NO new feature
generation or crosswalk needed — all on disk. Cohort-stratified Cox.

RESULT — FUSION_NEOANTIGEN -> OS: HR=1.24 per SD (95% CI 1.05-1.45, p=0.011),
i.e. MORE fusion-derived neoantigens -> WORSE overall survival. KM by
within-cohort fusion tertile: log-rank p=0.003, clean dose-response
(low>mid>high survival). SPLICE and ERV neoantigens: null for both OS and PFS.

RIGOUR CHECKS (signature-rigour-harness):
- Multiple testing: 6 tests, BH-FDR(fusion-OS)=0.066 — borderline, not <0.05.
- Cross-cohort replication: SAME DIRECTION in 3/4 cohorts (Gide HR1.54 p=0.023,
  Hugo 1.80 p=0.053, Liu 1.24 p=0.060; only Riaz reverses 0.75 n.s.). This is the
  strong part — not a one-cohort artifact.
- TMB confound: fusion is UNcorrelated with TMB (r=0.006). +TMB attenuates
  p 0.011->0.104 BUT univariate on the SAME n=191 TMB-subset is already p=0.094 —
  so the attenuation is POWER LOSS (72 samples, mostly Liu, lack TMB), NOT TMB
  confounding. HR barely moves (1.24->1.16). Signal is not a TMB proxy.

INTERPRETATION / CAVEATS (important, do not overclaim):
- DIRECTION IS OPPOSITE the antigenicity hypothesis. The project predicts more
  non-ref antigens -> more immunogenic -> BETTER ICB outcome. Fusion load goes the
  OTHER way (worse OS). So this is NOT confirmation of the evolutionary-RNA-state
  hypothesis as framed; if anything fusion-neoantigen load is a marker of a more
  aggressive/genomically-unstable tumor.
- PROGNOSTIC vs PREDICTIVE: these are single-arm ICB cohorts (no chemo/observation
  control), so this is an ON-TREATMENT survival association. Cannot distinguish
  "fusion load = generally worse prognosis" from "fusion load = worse ICB response
  specifically." A control arm (not available here) would be needed.
- These are DNA/WES-derived fusion neoantigen counts (a sibling pipeline), NOT the
  RNA fusion-transcript layer this campaign tried and failed to build on arm64.
  The RNA fusion layer remains uncomputed.
- Fusion is the ONE non-ref category with a survival signal; splice/ERV are null.

Artifacts: fusion_survival.png (8e00652b, v2 title qualified to "3 of 4 cohorts"),
fusion_os_forest.csv (ce32e10b).

### 9a. Fusion vs BINARY response/clinical-benefit (DONE — Mann-Whitney, per cohort)
Consistent with survival but underpowered vs the time-to-event test. Direction is
the SAME everywhere: in 3/4 cohorts responders carry LOWER median fusion load than
non-responders (median 1 vs 2). Per-cohort p: response gide 0.18 / hugo 0.084 /
liu 0.14 / riaz 0.12; clinical-benefit hugo 0.084 / liu 0.067 / gide,riaz n.s.
None individually significant, but the uniform direction matches the OS finding —
dichotomizing to RECIST just discards the timing info the Cox uses.

### 9b. MULTIVARIATE OS models — does fusion add value BEYOND immune / genomic load?
Key constraint: the joins are cohort-DISJOINT. Immune-floor scores are recoverable
ONLY for Gide (via data/registry/gide_id_crosswalk.csv, iatlas_sampleId<-run_accession,
n=59 PRE); TMB is available for everyone EXCEPT Gide. So no single model can hold
fusion + floor + TMB across cohorts. Two complementary models:
- MODEL 1 (Gide, n=59, 23 events): fusion + immune floor (Teff/TGFb) -> OS.
  Immune floor DOMINATES: HR 0.57 p=0.010 (more immune -> better OS, significant).
  Fusion attenuates 1.36->1.28 (p=0.25). BUT fusion was already n.s. univariate at
  this reduced n=59 (HR 1.36 p=0.14), so this does NOT cleanly show the floor
  "explains away" fusion — it is underpowered. Both keep their expected directions.
- MODEL 2 (non-Gide TMB cohorts liu+riaz+hugo, n=191, 100 events, cohort-stratified):
  fusion + TMB -> OS: fusion HR 1.16 p=0.104 (TMB HR 0.86 p=0.16). + SNV_neoantigen:
  fusion HR 1.15 p=0.15. Fusion barely moves under genomic-load adjustment and is
  UNcorrelated with TMB (r=0.006) -> NOT a mutational-load proxy.

VERDICT on independence: fusion's prognostic association is real univariate and
cross-cohort-directional, but MODEST and NOT cleanly shown to be INDEPENDENT of
immune composition — because in the only cohort where immune adjustment is possible
(Gide) the immune floor is the significant term and fusion attenuates (underpowered),
and where genomic-load adjustment is possible fusion survives only at p~0.10-0.15.
THE MISSING PIECE: immune-floor scores for Liu/Riaz/Hugo (needs scoring those
cohorts' expression) would allow the pooled fusion+floor model that actually settles
independence. That is the single highest-value next computation.

NEXT: (a) score the immune floor on Liu/Riaz/Hugo expression -> pooled fusion+floor
OS model (settles independence); (b) if the RNA fusion-transcript layer is ever
built, check concordance with this DNA fusion count.


## 10. DECISIVE pooled model (2026-07-12): fusion is INDEPENDENT of the immune axis
Ran the pooled fusion + immune-floor OS model — the composition-check methodology
applied to survival, and the test that resolves section 9b's open independence
question. Per sibling-agent methodological note: immune floor scored UNIFORMLY
across all 4 cohorts (same Ayers GEP-18 / IFNG-6 / Mariathasan Teff-10 / TGFb-23
gene sets, log2(x+1) -> within-cohort gene-z -> mean; salmon ENSG mapped to symbol
via GENCODE v46 GTF, Liu native symbols). All 18 GEP genes present in every cohort.
Floor scores saved /tmp/floor_uniform.parquet (259 samples).

JOIN: floor keyed run_accession(salmon)/Liu_SampleN -> analysis_frame sampleId via
gide_id_crosswalk (Gide), run_catalog sample_title regex Pt\d+ (Riaz +_pre / Hugo),
direct (Liu). POOLED FRAME n=252, 120 OS events, all 4 cohorts
(liu122/gide71/riaz32/hugo27).

POSITIVE CONTROL (validity gate): uniform floor predicts OS protective, stratified,
n=252: GEP HR 0.80 p=0.016; IFNG HR 0.80 p=0.014; Teff/TGFb HR 0.81 p=0.025. The
known-real predictor behaves correctly on this frame -> scoring + join are valid.

DECISIVE RESULT (cohort-stratified Cox, n=252, 120 ev):
  fusion + GEP:        fusion HR 1.288 (1.09-1.52) p=0.003 ; GEP HR 0.786 p=0.012
  fusion + Teff/TGFb:  fusion HR 1.269 (1.07-1.50) p=0.005 ; TGFb-bal HR 0.816 p=0.037
  fusion univariate (same n): HR 1.280 p=0.0036
=> FUSION AND THE IMMUNE FLOOR ARE INDEPENDENT, OPPOSITE-DIRECTION OS SIGNALS.
   Fusion HR is UNCHANGED (1.28->1.29) when the floor enters -> fusion is NOT immune
   composition in disguise, and (section 9b) NOT a TMB proxy. Two orthogonal signals:
   immune-protective, fusion-detrimental. This is the strongest finding of the campaign.

CAVEATS (unchanged, do not overclaim):
- Direction is OPPOSITE the antigenicity hypothesis (more fusion -> WORSE OS).
- Single-arm ICB cohorts -> prognostic (on-treatment), cannot prove ICB-predictive
  vs generally-prognostic without a control arm (not available).
- DNA/WES fusion-neoantigen COUNT, not the RNA fusion-transcript layer. Cross-layer
  replication (RNA fusions from the sibling STAR+Arriba run) is the next real advance
  and is NOT more proxy mining — per sibling coordination note, STOP further DNA-proxy
  slicing after this; wait for the RNA fusion output.
Artifacts: fusion_floor_multivariate.png (6d4b6e7d), .csv (5558cf3f), floor_uniform (/tmp).


## 11. RNA variant calling assessment + DNA variant-burden axes (2026-07-12)
Prompted by user Q on calling germline/somatic variants from RNA-seq + the Livne &
Yizhak FFixR paper (Bioinformatics 2026, RNA-MuTect-WMN + XGBoost FFPE-artifact filter,
validated on OUR cohorts Van Allen/Gide/Liu).

KEY FACTS:
- FFixR/RNA-MuTect-WMN is feasible + cohort-validated, BUT its headline validated
  output is RNA-TMB that RECOVERS DNA-TMB (R2=0.88). We already have DNA-TMB.
- IMPORTANT CAMPAIGN FACT: Gide + Liu are FFPE cohorts (per this paper). FFPE C>T/G>A
  deamination artifacts outnumber true mutations ~100:1 in naive RNA calls. A-to-I
  editing (A>G/T>C) is chemically distinct from FFPE C>T, so the editing layer is
  separable — but any RNA SNV work here must FFPE-filter.
- RNA-DNA discordance is large even frozen (>80% of RNA muts absent from DNA; only
  3.7-35% of DNA muts seen in RNA) -> RNA variant calls are intrinsically noisy.

DNA variant-burden -> OS (cohort-stratified, on-disk, NEVER tested before this):
  TMB_NONSYNONYMOUS  HR 0.775 p=0.0022   (protective — classic TMB->ICB)
  MUTATION_COUNT     HR 0.803 p=0.0052
  MUTS_CLONAL        HR 0.692 p=0.0019   (STRONGEST — truncal-neoantigen literature)
  MUTS_SUBCLONAL     HR 0.829 p=0.12  (ns)
  SNV_NEOANTIGEN     HR 0.932 p=0.53  (ns) ; INDEL_NEOANTIGEN 0.896 p=0.32 (ns)
  -> mutation/clonal BURDEN is protective; per-category SNV/indel neoantigen counts null.

3-AXIS MODEL fusion + TMB + immune floor -> OS (n=180, no Gide [Gide lacks TMB], 91 ev):
  fusion HR 1.218 p=0.043 (survives) ; tmb HR 0.90 p=0.37 ; imm HR 0.90 p=0.38 ; C=0.618
  INTERPRETATION (do not overclaim): fusion is the most robust — survives immune
  adjustment (n=252), TMB adjustment, and the joint model. TMB + immune attenuate here
  but this frame DROPS Gide (immune's strongest cohort) and is Liu-heavy, so their
  attenuation is partly power/subset, NOT demonstrated redundancy. Both are significant
  univariate at full n (TMB n=335 p=0.002; immune n=252 p=0.012).

RECOMMENDATION ON RNA VARIANT CALLING: NOT worth the heavy RNA-MuTect-WMN/FFixR build
for this project. Its validated value (RNA-TMB) recovers DNA-TMB we already have and
which already works. The only non-redundant angles are (a) EXPRESSED-mutation TMB
(Katzir reports modest gains) and (b) RNA-DNA discordance as an immunoediting readout
(DNA mutation not expressed = possible immune escape). Both are targeted refinements,
not classifier features, and both need the FFPE filter. Reserve for if the fusion
finding needs mechanistic follow-up.
Net: campaign now has 3 prognostic axes — immune (protective), mutation/clonal burden
(protective), fusion neoantigen (DETRIMENTAL, most robustly independent). Continue to
hold DNA-proxy mining per section 10; RNA fusion cross-layer replication remains the
next real advance.


## 12. Expressed-TMB test (2026-07-12): the RNA-variant angle, resolved empirically
Ran the one genuinely non-redundant RNA-variant angle the FFixR paper names —
EXPRESSED-mutation TMB (Katzir/Sorokin claim RNA-based TMB reflecting only expressed
mutations has "improved predictive power"). No RNA variant CALLING needed: fetched the
full per-mutation MAF for the 3 iAtlas cohorts with MUTATION_EXTENDED (Liu 78k, Riaz
20k, Hugo 39k rows, cBioPortal REST /mutations/fetch), then for each sample counted
mutations falling in genes EXPRESSED in that sample's own RNA (TPM>1; Liu cbioportal
symbol expr, Riaz/Hugo salmon ENSG->symbol via GENCODE v46). Gide excluded (no iAtlas
mutation profile — consistent with its missing frame-TMB).

VALIDITY: MAF-derived nonsyn DNA-TMB vs analysis_frame TMB_NONSYNONYMOUS r=0.992 (n=179)
-> pipeline correct.

RESULT (cohort-stratified Cox, n=179, 91 OS events):
  total DNA-TMB   HR 0.903 (0.73-1.12) p=0.36  C=0.559
  expressed-TMB   HR 0.914 (0.73-1.14) p=0.43  C=0.557
  Liu-only (n=122): dna 0.816 p=0.12 vs expr 0.820 p=0.14 — identical
=> EXPRESSED-TMB IS PROGNOSTICALLY INDISTINGUISHABLE FROM TOTAL DNA-TMB.
WHY: 77% of DNA mutations fall in expressed genes and the two measures correlate r=0.98
(near-collinear). Restricting to expressed mutations removes little and changes nothing.

CONCLUSION: this closes the RNA-variant question empirically, not just by argument. The
expressed-mutation refinement — the single non-redundant thing RNA variant calling could
have added on THIS data — provides no prognostic gain over the DNA-TMB already on disk.
Combined with the FFPE-artifact hazard (Gide/Liu FFPE) and heavy RNA-MuTect-WMN/FFixR
build cost, RNA somatic variant calling is NOT a worthwhile direction for this project.
(Note this contradicts nothing in the paper — the paper's win is recovering DNA-TMB when
DNA is UNAVAILABLE; we have DNA, so the recovery has no marginal value here.)
The other RNA-variant angle (RNA-DNA discordance as immunoediting readout) genuinely
needs RNA variant calls and is deferred; it is a mechanistic follow-up, not a classifier
feature. Standing conclusion unchanged: 3 prognostic axes (immune protective / mutation
burden protective / fusion detrimental+robust), hold DNA-proxy mining, RNA fusion layer
is the next real advance.
Artifacts: expressed_tmb.png (9c9d33ab), expressed_tmb.csv, /tmp/expressed_tmb_os.parquet.


## 13. Incremental value OVER TMB — the decision-relevant test (2026-07-12)
Q (user): in addition to TMB, do any annotations ENHANCE prediction? Nested Cox LRT on
a FIXED sample set (n=180, 91 OS events, liu122/riaz32/hugo26; all candidate features
non-null so the LRT is valid), cohort-stratified. Base = TMB; also base = TMB + immune floor.

RESULT (incremental over TMB alone):
  + FUSION_NEOANTIGEN : LRT p=0.058, dC +0.031, HR 1.21 p=0.053   -> borderline ENHANCES
  + SPLICE_NEOANTIGEN : LRT p=0.979, dC  0.000                     -> NO value
  + immune floor (gep): LRT p=0.507                                 -> no value over TMB here (subset)
Incremental over TMB + immune floor:
  + FUSION_NEOANTIGEN : LRT p=0.048, dC +0.041, HR 1.22 p=0.043   -> ENHANCES (significant)
  + SPLICE_NEOANTIGEN : LRT p=0.954                                 -> NO value
Trustworthy full model TMB + floor + fusion: fusion HR 1.218 p=0.043, C=0.618.

ANSWER: FUSION is the ONLY annotation that adds prognostic value beyond TMB (and beyond
TMB+floor). Splice-neoantigen count adds nothing. So the incremental story is:
mutation burden (TMB) + fusion — fusion is the one non-TMB feature that enhances.

*** ERV ARTIFACT — DO NOT RE-CHASE (rigour harness catch) ***
An initial run showed ERV_NEOANTIGEN "enhancing" over TMB (LRT p=0.023, HR 1.27-1.31
p=0.002-0.006), which looked like a TE-activation win supporting the hypothesis. IT IS
SPURIOUS. ERV_NEOANTIGEN is near-constant ZERO across the whole dataset: only 5/264
samples are nonzero (Gide 2 [value 1], Liu 3 [two value 1, one value 2]; Hugo/Riaz all 0).
In the n=180 modeling frame (liu/riaz/hugo — Gide excluded, no TMB) the 3 nonzero are all
Liu, so the cohort-stratified Cox was keying on THREE individual Liu patients. It also breaks the
within-cohort permutation (zero-variance strata). DISCARD ERV_NEOANTIGEN as a WES-derived
feature — the DNA annotation pipeline essentially never called ERV neoantigens. (This is a
DATA-SPARSITY artifact, NOT evidence about TE biology — TE/ERV must be tested via the
RNA family-level TE quantification, not this empty WES column.) ERV vs TMB r=0.03, ERV vs
fusion r=-0.02 (orthogonality was real, but on 3 nonzero points it is meaningless).
Splice-neoantigen, by contrast, is genuinely variable (59-64 distinct/cohort) and its null
is a REAL null.
Artifacts: incremental test in-kernel (/tmp/pooled_fusion_floor.parquet + analysis_frame).
