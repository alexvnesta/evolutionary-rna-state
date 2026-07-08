# Move 2 + hypothesis-kill tests — watch & auto-execute trigger

## Status: WIRED and RUNNING on partial data (68 samples, 3 cohorts)
quant_gene_tpm.parquet now holds **68 pre-treatment samples (gide2019: 30, hugo2016: 28,
riaz2017: 10)**, fully populated, 62,266 ENSG genes. hugo landed mid-run and the watch
auto-fired; the hugo run->iAtlas crosswalk was built and validated (26/28 map; ENA
PRJNA312948, sample_title IS the iAtlas patientId 'Pt##', with Pt27A/B -> Pt27; Pt16 absent
from the iAtlas response set). The pipeline is still filling cohorts (liu/dfci pending).
The runner RUNS now and re-fires as the matrix grows. **liu2019 (n=122, has TMB) is the
decisive remaining input** — it takes the response tests to n~180 and a 4th held-out cohort.

## The auto-runner now fires FOUR test families on every matrix change
move2_autorun.py v7 (run_if_changed) detect->score->crosswalk->test->save, plus:
  1. T1/T2/T3 regulator-activity tests (run_activity_response_tests).
  2. reframing P1-P3 (reframing_tests): coordination null, immune-cold coupling, IFN-incremental.
  3. **learned-representation head-to-head** (headtohead_repr): PCA embedding vs IFN immune
     floor under random-5fold AND leave-one-cohort-out, with a LOCO permutation null. This is
     the untested-core experiment — the first honest test of the hypothesis's stated method.
  4. per-cohort ID-mapping diagnostics that FLAG any cohort whose IDs don't resolve to iAtlas
     (rather than silently dropping them) — this is how the hugo/liu ID mismatches were caught.

## What was built this session
- move2_autorun.py — idempotent detect->score->crosswalk->test->save runner.
  * Detects latest quant_gene_tpm.parquet; short-circuits (11 ms) if version_id unchanged.
  * Scores regulator activity (rbp_activity_scorer, 56 genes -> 3 sets), all 56 ENSG
    present in the matrix; HGNC->ENSG via Ensembl REST (cached to .regulator_ensembl_map.json).
  * Builds per-cohort run->iAtlas crosswalk: gide from the validated artifact; riaz derived
    from ENA PRJNA356761 sample titles (Pt##_Pre_... -> iAtlas Pt##). VALIDATION of the
    10 riaz samples in the current 40-sample matrix (the run->Pt## mapping): all 10 join
    1:1 (unique run, unique Pt), all 10 are PRE-treatment, all 10 resolve in the iAtlas
    riaz frame with RESPONSE populated, and all 10 are independently triangulated in the
    Riaz mmc3 mutation supplement (built from unrelated sample IDs) -- so the Pt## identity
    is corroborated by a second, independent source. NB: this is ID-identity validation;
    it does not include gide-style RECIST-vs-arm concordance (riaz titles carry Pt## and
    timepoint literally, with no arm-coding rule to concordance-check).
  * Runs run_activity_response_tests (T1 mechanistic / T2 shared-state / T3 LOCO-AUROC).
  * Writes results_move2_<n>samp_* and a .move2_state.json fingerprint.
- move2_watch_live.py — bounded in-session poll loop (fires the runner on version change).
- run_activity_response_tests.py — UPDATED (v2). Two correctness fixes:
  1. crosswalk columns (arm/timepoint/cohort) now OPTIONAL, not required.
  2. **join is now cohort-scoped** (cohort, iatlas_patientId). Critical: hugo & riaz
     share 'Pt##' patient IDs (14 collisions), so a cohort-blind join created phantom
     cross-cohort matches. Now fixed.

## HOW TO RUN (any future session, one call, in the erna kernel)
    import importlib.util
    spec = importlib.util.spec_from_file_location("m2","move2_autorun.py")
    M2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(M2)
    print(M2.run_if_changed(host))          # no-ops if matrix unchanged; runs if grown
    # then: save_artifacts(["results_move2_<n>samp_report.json", "..._T1_mechanistic.csv"])
Dependencies (all in project store, copy into workspace first if missing):
    rbp_activity_scorer.py, run_activity_response_tests.py (v2), reframing_tests.py,
    headtohead_repr.py, analysis_frame.parquet, gide2019_id_crosswalk.csv
Adding a NEW cohort's expression: gide/riaz/hugo have validated ENA-title rules in
build_crosswalk(). liu2019 arrives "via DFCI matrix (processed)" and is handled by the
_passthrough_crosswalk() generic path, which resolves matrix IDs against BOTH the iAtlas
patientId AND sampleId columns. THIS IS LOAD-BEARING for liu: its iAtlas records key on
patientId='Patient4' but the processed matrix carries sampleId='Liu_Sample4' — verified
122/122 liu samples map to the patientId join key under EITHER form. Any cohort whose IDs
match neither is emitted as-is and FLAGGED by the diagnostics (not silently dropped).

## LIMITATION on the watch ("auto-execute")
There is no cross-session daemon: a live watch only runs while THIS session's kernel is
alive (move2_watch_live launched here, 3-min poll, 6-h self-terminating cap). For a
durable trigger, a future session must call run_if_changed(host) once — it is idempotent
and cheap, so it is safe to drop at the top of any pipeline-checkpoint session.

## Current result (66 samples, gide+hugo+riaz) — the headline finding
- **Learned-representation head-to-head (the untested core)**: with a genuine 3rd held-out
  cohort, the learned RNA representation is at CHANCE across cohorts (LOCO AUROC 0.49, perm
  p=0.41), while the immune floor is real-but-weak (0.59, perm p=0.044). Combining them drags
  the floor to chance (0.50). Random-5fold parity (rep 0.75 vs floor 0.76 at n=40) was
  within-cohort overfitting. 9-config hyperparameter sweep holds the rep LOCO in 0.49-0.63.
  VERDICT: a learned RNA representation carries no cross-cohort ICB signal beyond immune
  composition and degrades the predictor that works. (headtohead_repr_66samp.json.)
- **Residual tumor-intrinsic over hardened immune floor (n=66)**: joint p=0.383 (nothing
  survives); immune floor over intrinsic p=0.031. Cohort-adjusted, 4/11 immune cell types
  significant (checkpoint/CD8/Treg/CD4 — the T-cell axis), stroma/endothelial flat.
- Activity-alone tests (T1/T2/T3) remain underpowered leads, not results.
INTERPRETATION: the hypothesis's antigen-quantity core is dead (416 samples, prior work), and
its learned-representation method — tested honestly for the first time — fails at cross-cohort
transfer. What survives is immune-composition dominance (the CIBERSORTx critique confirmed) and
a DNA-side heterogeneity effect. liu (n=122, has TMB) is the powered adjudication: 4th held-out
cohort, ~n180. The current n=40->66 trajectory strongly predicts the powered tests stay null.

## Artifacts (current)
- move2_autorun.py (v7), move2_watch_live.py, run_activity_response_tests.py (v2),
  headtohead_repr.py, reframing_tests.py
- combined_id_crosswalk.csv (v2, 3 cohorts), results_move2_68samp_report.json
- headtohead_repr_40samp.json, headtohead_repr_66samp.json,
  residual_over_hardened_floor_40samp.json, residual_over_hardened_floor_66samp.json,
  learned_rep_sensitivity_40samp.json, immune_deconvolution_block_40samp.csv
- FINAL_VERDICT.md, fig_master_verdict.png (top-level consolidated read)
