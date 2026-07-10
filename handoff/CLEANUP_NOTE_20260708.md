# Repo cleanup note — 2026-07-08

For sibling agents / future sessions. A hygiene pass was run on the repo root.
This records **what changed**, **what was deliberately left alone**, and **one
pending task** so nobody redoes it or gets surprised. Nothing here touches the
science; see `CLAUDE.md` for the live cross-agent working notes (unchanged).

## Done (disk reclaim — ~150 GB out of the working tree)

Deleted (moved to Trash) the Nextflow **work scratch** for completed, published
runs under `results/large/`:

- `nf_work_rnaseq` (149 GB), `nf_work_ir`, `nf_work_editing{,2,3}`.
- `nf_work_preview` (0 B, empty directory) — removed as empty clutter, not
  covered by the output-verification below (nothing to lose).

Safe because published outputs were verified intact **outside** scratch first:
- `results/rnaseq_pilot_hisat2/` — full nf-core output; `PD1_35_PRE.markdup.sorted.bam` BGZF-complete (not truncated).
- `results/editing_bams/` — all 16 HISAT2 BAMs BGZF-complete (the real editing-arm data; **not** scratch, untouched).
- Only `-resume` caching for those specific runs is lost. Re-runs start fresh.

Also removed transient git-ignored clutter at root: rotated `.nextflow.log.{1..9}`
and `__pycache__/`. The current `.nextflow.log` was kept.

> Space note: deletes go to Trash on the same volume — **empty the Trash** to
> actually reclaim the ~150 GB.

## Left ALONE on purpose

- **`.nextflow_home/` (10 GB) — KEEP.** It is the staged *offline* execution
  environment (pre-downloaded plugins because the registry is proxy-blocked,
  GNU-tools env, local conda channel, persistent tximeta env), not scratch.
  Deleting it forces rebuilds through blocked registries. See `pipelines/RUNBOOK.md`.
- **`results/large/nf_work_te_erv` — KEEP.** A TE/ERV (Telescope) run was in
  flight (log showed `TELESCOPE_ASSIGN` running). This is active/high-value
  work, not a completed run — do not clear it until it is confirmed done and
  its outputs are published to `results/te_erv/`.
- **`CLAUDE.md` — not edited.** It is the live cross-agent doc owned by the
  active branch session; updates to it belong to that session.

## PENDING — root shell scripts should move to `pipelines/scripts/` (not done yet)

Seven `.sh` helpers sit loose at repo root, violating the project convention
that shell helpers live in `pipelines/scripts/`. The move was **deferred** to
avoid colliding with active work: `align_hisat2.sh` and `run_aei_panel.sh` have
**uncommitted edits** from the branch session (BAM → lossless-CRAM migration in
progress), and `archive_bam_to_cram.sh` was running. Moving/committing files
another agent is mid-edit on would cause exactly the confusion we want to avoid.

Proposed relocation once the CRAM migration is committed (all use absolute
`$REPO/...` paths internally, so moving them does not change behaviour; only the
git tracking + a few doc references need updating):

| Script | Action | Why |
|--------|--------|-----|
| `align_hisat2.sh` | → `pipelines/scripts/` | current aligner (has uncommitted CRAM edits) |
| `run_aei_panel.sh` | → `pipelines/scripts/` | current AEI driver (has uncommitted CRAM edits) |
| `prefetch_fastqs.sh` | → `pipelines/scripts/` | current, referenced in docs |
| `align_editing_subset.sh` | → `pipelines/scripts/` | current subset helper |
| `align_only.sh` | retire (or `pipelines/scripts/legacy/`) | STAR-based; STAR is broken on arm64 (RUNBOOK) |
| `run_aei_batch.sh` | retire (or legacy/) | superseded BAM-era AEI |
| `run_aei_parallel.sh` | retire (or legacy/) | superseded BAM-era AEI |

Doc references to fix after the move: `docs/PROJECT_STATUS.md` (lines ~74–75),
`docs/PIPELINE_HANDOFF.md` (lines ~87, ~104). `compute_aei_fast.py` is called by
`run_aei_panel.sh` via `$REPO/compute_aei_fast.py` (absolute), so it can move or
stay independently — but keeping the script and the python it calls together in
`pipelines/scripts/` is tidier.

Whoever picks this up: git dir is `.gitmeta` — drive with
`GIT_DIR=.gitmeta GIT_WORK_TREE=. git mv <src> pipelines/scripts/`.

## Hygiene scan (read-only, `evolutionary-rna-state-cleanliness` skill)

Score 94/100. Open items unrelated to the above: `hugo2016_S1D.csv` (11.6 MB
committed — consider LFS/ignore), `analysis/clonal_trajectory/code` lacks a
README, and two untracked source items (`pipelines/scripts/archive_bam_to_cram.sh`,
`pipelines/te_erv/samplesheets/`) — commit once the branch session's work settles.
