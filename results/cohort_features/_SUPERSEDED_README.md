# results/cohort_features/ — SUPERSEDED (2026-07-10)

This directory was a SECOND, redundant non-reference feature build (session 756ba1cd,
`cohort_driver.sh`), computing the same TE-family / intron-retention / splicing-junction /
AEI features as the canonical build in `results/nonref_run/` (session 15defe54).

**Consolidation decision (redundancy audit 2026-07-10):**
`results/nonref_run/` is the SINGLE canonical non-reference feature build. It aligns all
n=106 samples with a uniform `hisat2 -k 10` (multimappers retained) — the correct substrate
for TE/repeat quantification. This cohort_features build mixed default-k reused BAMs with
plain-HISAT2 downloads, which would put an internal alignment batch effect into the feature
matrix that gets tested for ICB signal.

The driver was stopped gracefully (per-sample STATUS skip-markers -> loop reached EOF/COMPLETE
at 23:24 on 2026-07-10). The 5 genuinely-built samples here (ERR2208909, SRR3184285,
SRR3184288, SRR5088840, SRR5088867, SRR5088891) are LEFT INTACT for cross-checking against the
canonical run but are NOT the project's feature source of record.

Do not relaunch cohort_driver.sh. To extend the cohort, drive results/nonref_run/ instead.
See docs/PROJECT_STATUS.md and handoff/AGENT_REDUNDANCY_AUDIT_20260710.md.
