# Active-agent redundancy & hygiene audit — 2026-07-10 02:45 UTC
Auditor session: 64079601 (user-launched). Scope: all `processing` frames in project proj_98c1e1f16a1a.

## 5 sessions currently PROCESSING

| Session | Nominal role | What it is ACTUALLY doing now | Evidence |
|---|---|---|---|
| 837512d2 Raw RNA-seq Feature Fine-tuning | (repurposed) cross-session **coordinator** | Owns `docs/PROJECT_STATUS.md`; polls all frames, commits status snapshots, monitors pipelines | 451 cells; last cells all `host.query`/`host.frames` polling + snapshot writes |
| 756ba1cd Project Status Overview & Prioritization | forensic synthesis + prioritization | (1) Unique 4-audit forensic re-audit (hypothesis UNTESTED not falsified; EVA=expr-PCA R²=1.0; Orthrus never run at scale; mediation analysis). (2) NOW ALSO running the non-ref caller pipeline in `results/cohort_work/` | 334 cells; cohort driver launched 22:39 |
| 64079601 Audit Active Agents (THIS) | redundancy + doc/commit audit | This report | 10 cells |
| c15a540e Find/Create Nextflow Pipelines | raw-read pipeline (arm64) | nf-core/rnaseq **HISAT2 alignment** of Gide cohort FASTQ→BAM | 975 cells (most active); 22 `*.markdup.sorted.bam` in `results/rnaseq_cohort/hisat2/` |
| 15defe54 Evaluate Orthrus | non-ref feature builder | `results/nonref_run/`: own HISAT2 align + AEI/TE-family/IR/junction callers | 366 cells; manifest_n106 (all existing_bam=0) |

## REDUNDANCIES

### A. Non-reference callers run twice on the same sample (DEFINITE, live)
- Orthrus (15defe54) COMPLETED `SRR5088840` at 18:40 → `results/nonref_run/out/SRR5088840/` has te_family/ir_wide/junctions/aei + STATUS.
- Overview (756ba1cd) RE-STARTED the identical TE-family + junctions + IR callers on the SAME `SRR5088840` at 22:39 → `results/cohort_work/SRR5088840/`, currently churning ~7.9 GB of featureCounts `temp-core-*.sam` scratch.
- Two sessions writing overlapping non-ref outputs to different dirs → wasted CPU + risk of divergent feature matrices.

### B. FASTQ→BAM alignment duplicated for the Gide cohort (LIKELY)
- Nextflow (c15a540e): nf-core/rnaseq HISAT2 → `results/rnaseq_cohort/hisat2/*.markdup.sorted.bam` (trial-ID names PD1_*/ipiPD1_* = Gide).
- Orthrus nonref_run manifest_n106: all 106 rows `existing_bam=0` → aligns every sample itself with its own HISAT2, incl. the same Gide ERR22089xx runs.
- Two independent HISAT2 passes over the same Gide FASTQs. No shared BAM contract between the two pipelines.

### C. Coordination/status tracking triplicated
- 837512d2 (owns PROJECT_STATUS.md), 756ba1cd (forensic synthesis/prioritization), 64079601 (this audit) all consolidate project state. PROJECT_STATUS.md itself flags "Audit Active Agents ... overlaps this doc."

## DOCUMENTATION
- Strong overall: PROJECT_STATUS.md is a thorough live coordination doc; each session self-commits; forensic corrections captured.
- GAP: new BCR/SHM subsystem (`analysis/differentiated/bcr_shm.{py,nf}`, `bin/`, `modules/`, `test_bcr_shm.py`, `pipelines/bcr_repertoire/run_trust4_pilot.sh`) — substantial Jul-9 work implementing the c8f0a686 SHM advisory — is NOT referenced in any committed doc (README/CLAUDE/docs).

## COMMITS  (GIT_DIR=.gitmeta; HEAD == origin/main)
26 uncommitted entries:
- Modified (4): analysis/pilot_gep.py, docs/ENCODER_REVIEW.md, docs/FINAL_VERDICT.md, handoff/HUGO_COMPLETION_28.md
- Untracked (22): BCR/SHM subsystem; analysis/baseline/{apm_scores,tls_bcell_scores}(+tests); analysis/pilot/run_salmon_*.sh; data/manifests/selection_manifest_remaining.csv; handoff/CLEANUP_NOTE_20260708.md; artsave/; tools/; CLAUDE.md
- tools/ holds large third-party binaries (STAR-2.7.11b(.tar.gz), TRUST4-1.1.5, arcasHLA) — should be .gitignored, NOT committed.
- CLAUDE.md and claude.md: BOTH untracked, byte-identical, same mtime — case-collision on the case-insensitive macOS FS (git ignorecase=true). Latent hazard.
