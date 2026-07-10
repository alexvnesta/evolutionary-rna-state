# Canonical non-reference feature pipeline — corrected decision (2026-07-10)

_Author: audit/coordination session 64079601. This SUPERSEDES the 2026-07-10 ~23:25
consolidation banner in PROJECT_STATUS.md, which picked the wrong survivor._

## The correction

An earlier consolidation crowned the ad-hoc bash build `results/nonref_run/`
(session 15defe54) as the single canonical non-reference feature build, on the
stated reasoning that its uniform `hisat2 -k 10` alignment was "the correct
substrate" and the Nextflow pipeline's per-feature alignment was a "mixed regime
→ batch effect." **That reasoning was inverted and the decision was wrong.**

Per-feature alignment is the *correct* methodology, not a defect:
- **RNA editing (AEI/JACUSA2)** needs UNIQUE mappers (MAPQ 60) — multimappers
  create false A-to-I calls.
- **TE/ERV locus (Telescope)** needs MULTIMAPPER-permissive alignment
  (`bowtie2 -k 100 --very-sensitive-local`, Bendall et al. 2019) — reassigning
  ambiguous reads to specific loci is Telescope's whole purpose.
- A single `-k 10` base alignment for all features is a COMPROMISE, not a gold
  standard.

## Decisive evidence: the bash build is strictly less complete

`results/nonref_run/` **dropped locus-level Telescope entirely.** Its own driver
notes (and PROJECT_STATUS.md §298-299) record that Telescope on the HISAT2 `-k10`
BAM failed twice (`KeyError: 'locus'`, then `TypeError` in `_load_sequential`),
the bowtie2 retry was interrupted at ~95 min before Telescope ran, and locus-TE
was then dropped in favor of **family-level TE via featureCounts** only.

The comprehensive Nextflow pipeline (`pipelines/`, session c15a540e) is
**designed** to run locus-level Telescope on a `bowtie2 -k 100` multimap pass
with a patched arm64 Telescope v1.0.3 wheel — a feature layer the bash build
explicitly dropped.

**IMPORTANT — verification caveat (checked on disk 2026-07-10):** the RUNBOOK
(`pipelines/RUNBOOK.md §3`) narrates a completed te_erv integration run (exit 0)
producing a 15,209-locus matrix + 1,328-family matrix at
`results/te_erv_integration_test/matrices/`. **That output could NOT be verified
on disk in this session:** the directory `results/te_erv_integration_test/` does
not exist, no `*telescope_report.tsv` / `*.cntTable` / locus matrix files exist
anywhere under `results/`, the Nextflow work-dirs for TELESCOPE_ASSIGN and
MERGE_TELESCOPE were not found (cleaned or never completed), and the only te_erv
artifact on disk is a 13-line `results/te_erv_integration_test.log` in which
TELESCOPE_ASSIGN and MERGE_TELESCOPE are only "Submitted" — with no completion
banner. So locus-level TE is a **documented-but-unverified** capability: the code
and the patched wheel are in the repo, but I have NOT confirmed a successful
locus-TE run produced real output. Treat the RUNBOOK's te_erv numbers as an
unverified prior claim, not established fact, until a cohort run reproduces them.

What IS verified on disk for TE: **family-level** counts only
(`results/phase0_proof/ERR2208909/te_family_counts.txt`, and the bash build's
`te_family.counts`). Both pipelines therefore have *verified* family-level TE;
**neither has verified locus-level TE output on disk** — the Nextflow pipeline
has the code path for it and the bash build removed it.

## Verified pipeline status (RUNBOOK §3 + on-disk output census, 2026-07-10)

| Phenotype | Pipeline | Real output verified? |
|---|---|---|
| Base align + quant | nf-core/rnaseq 3.26.0 (HISAT2+Salmon, arm64) | YES — spine BAMs + Salmon |
| Alt. splicing | nf-core/rnasplice 1.0.4 (DEXSeq+edgeR; genome_bam) | YES — pilot exit 0, DEXSeq 74 sig bins |
| RNA editing | custom `pipelines/rna_editing/` (JACUSA2 + AEI) | YES — `results/editing_integration_test/` aei tsv |
| Intron retention | custom `pipelines/intron_retention/` (featureCounts IR) | YES — `results/ir_integration_test/` ir tsv |
| TE / ERV | custom `pipelines/te_erv/` (Telescope locus + TEcount family) | FAMILY-level: YES on disk (`phase0_proof/.../te_family_counts.txt`). LOCUS-level (Telescope): code present, **NOT verified on disk** — RUNBOOK claims a run but no matrices/work-dirs exist (see caveat above) |
| Fusion transcripts | nf-core/rnafusion 4.1.3 (Arriba + STAR-Fusion) | NO — arm64-BLOCKED (needs STAR; runs on amd64/Docker only) |
| Composition / purity | InstaPrism / NMF deconvolution (staged `tools/r-deconv-libs`) | Ran in analysis/ (kill_tests immune block, clonal DDLPS); NOT yet wired as an NF module |

Note on the te_erv evidence trail (to prevent future over-claims in EITHER
direction): this session first said te_erv had "only a log, no output," then
over-corrected to assert the RUNBOOK's completed locus/family matrices as fact.
Both were wrong. The disk truth (verified 2026-07-10): family-level TE output
exists (`phase0_proof/.../te_family_counts.txt`); locus-level TE (Telescope)
output does NOT exist on disk and its integration run could not be verified —
only a 13-line log with "Submitted" steps. The code path and patched wheel are
present; the run is unproven here. Fusion is separately arm64-blocked.

## DECISION

**Canonical = the comprehensive arm64 Nextflow pipeline** (`pipelines/`, session
c15a540e): nf-core/rnaseq spine → {rnasplice, rna_editing, intron_retention,
te_erv} fan-out, all on the HISAT2 arm64 spine, per RUNBOOK.md.

This decision rests on grounds that ARE verified on disk this session, plus the
user's explicit direction:
1. **Correct per-feature alignment (verified from both pipelines' scripts):**
   editing uses unique-mappers (MAPQ 60), TE-locus is architected for a
   multimap `bowtie2 -k 100` pass — vs the bash build's single `-k 10` regime
   for everything. This is a methodology fact independent of any run.
2. **Verified-on-disk feature output:** rnaseq spine BAMs+Salmon, rnasplice
   (DEXSeq exit 0), rna_editing (aei tsv), intron_retention (ir tsv),
   family-level TE. These are real outputs.
3. **User directive (2026-07-10):** "Nextflow pipeline canonical, stop bash
   nonref_run."

It is NOT (yet) fully proven end-to-end: locus-level TE is unverified on disk
(see caveat), and cohort-scale runs of the fan-out arms remain to be done. The
remaining work is cohort scale-up PLUS a real locus-TE verification — not
greenfield pipeline development.

**Deprecated = `results/nonref_run/`** (bash compromise). Stopped 2026-07-10 (see
PROJECT_STATUS.md). Its completed per-sample outputs are kept for cross-checking
only; they are NOT the feature source of record (no locus-TE, single-regime
alignment).

**Deferred = fusion + deconvolution-as-NF-module.** Fusion needs an amd64/Docker
host. Deconvolution has run in analysis/ but is not yet a pipeline module; fold
it in as an explicit composition/purity step (it is named in ROADMAP.md).
