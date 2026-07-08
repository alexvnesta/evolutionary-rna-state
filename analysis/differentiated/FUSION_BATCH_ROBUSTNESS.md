# Fusion neoantigen burden — batch robustness

FUSION_NEOANTIGEN_BURDEN — batch/platform robustness
====================================================
The MHC-binding half of this feature is batch-invariant (MHCflurry percentile
rank is calibrated per allele against a fixed random-peptide background, so it
does not drift with input composition). The FUSION-DETECTION half is NOT: fusion
call sets depend strongly on the caller, its version, the reference/annotation,
and — critically — sequencing depth and read length (more reads -> more junction
support -> more low-confidence calls). Uncontrolled, this injects a technical
axis that mimics biological signal.

FIXED-CALLER REQUIREMENT (enforced in code):
  * `fusion_features_for_sample` accepts exactly ONE caller input (Arriba XOR
    STAR-Fusion) and records `caller` + `caller_version` on every row. Never
    mix Arriba-derived and STAR-Fusion-derived burdens on the same comparison
    axis; the two callers have different sensitivity/precision profiles and
    different in-frame annotation logic.
  * Hold the caller VERSION and the reference/annotation build constant across
    every sample in a comparison (record both).

CONFIDENCE FILTER: default to Arriba high+medium confidence (drops the depth-
sensitive low-confidence tail). This is the main lever against depth-driven
false positives. STAR-Fusion analogue: require FFPM / junction+spanning read
support thresholds and `--examine_coding_effect`.

DEPTH: where library sizes differ substantially across a cohort, either
down-sample to a common depth before calling, or include per-sample fusion read
depth as a covariate. Report `n_fusions` and `n_inframe_fusions` alongside the
burden so a depth confound is visible (a cohort with systematically higher
n_fusions likely differs in depth, not biology).

REPORTING: z-score / rank within cohort before pooling across cohorts, exactly
as for the sibling burdens; never pool raw counts across platforms.
