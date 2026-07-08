# Data-inventory → pipeline handoff: cohorts, joins, and the quantification queue

Status brief from the data-inventory session for whoever runs the RNA
quantification pipeline. Everything below is committed to the repo and
verified against source; identifiers and counts are copied from the saved
files, not from memory. Nothing here is a plan to re-derive. The data
foundation is complete, and the remaining work is compute (quantify the full
run set), not data assembly.

## The binding limitation

**The 3-cohort transfer test currently rests on n=10–12 pilots, and that is the
binding limitation on every downstream claim.** The Gide-trained de-novo antigen
axis fails to transfer to both held-out cohorts (AUROC 0.36–0.44 on Riaz n=10,
0.56–0.58 on Hugo n=12 with 95% CI [0.22, 0.92] spanning chance), versus ~0.87
within Gide. But the Hugo and Riaz arms are pilots, so "non-replication" cannot
yet be separated from "underpowered." Source: `results/heldout_transfer_3cohort.json`
and `results/deepen_results.json`. The single highest-value pipeline action is to
quantify the full run set below so those held-out arms leave the pilot regime.

## What is cataloged and ready to quantify

`data/catalog/run_catalog.csv` — **228 RNA-Seq runs**, all three RNA cohorts,
with FASTQ FTP URLs + MD5s and joined clinical fields (RECIST, responder/
non-responder, therapy, timepoint, OS, vital). Quantification queue by cohort:

| Cohort   | Runs | Data volume | Pilot done? | Notes |
|----------|------|-------------|-------------|-------|
| riaz2017 | 109  | ~486 GB     | partial     | Illumina Genome Analyzer; more read pairs, ~half the bases of Gide |
| gide2019 | 91   | ~824 GB     | pilot (n≈30 used in transfer) | HiSeq 2500; deepest/longest reads; RNA development lead |
| hugo2016 | 28   | ~239 GB     | **only n=12 quantified** (`results/hugo_gene_tpm.parquet`) | HiSeq 2000; 16 runs still unquantified |

**Immediate target: the 16 un-quantified Hugo runs.** `hugo_gene_tpm.parquet`
currently holds 12 samples; the catalog now defines all 28 (SRA study SRP070710
/ BioProject PRJNA312948, GSE78220, verified human RNA-Seq via ENA GEO→SRA
elink). Completing Hugo moves the weakest arm of the transfer test (95% CI
[0.22, 0.92]) toward a decidable result first, at the lowest data cost of the
three cohorts.

## Joins that are built, validated, and waiting on quantified expression

These unblock the two analysis directions named in `HANDOFF_rna_state_next.md`.
Full method + caveats: `data/registry/_joins_provenance.json`.

1. **Gide expression ↔ iAtlas burden/neoantigen features** —
   `data/registry/gide_id_crosswalk.csv`. Run-level map from ENA PRJEB23709 IDs
   (`ipiPD1_N` / `PD1_N`) to iAtlas `sampleId` (`iPiN_On` / `PD0N_Pre`).
   91 runs / 75 patients, validated by clinical fingerprint (75/75 concordant on
   arm, timepoint structure, and response; 91/91 runs round-trip to valid iAtlas
   IDs). Lets the RBP/splicing-factor-activity scorer link expression to the
   SPLICE/ERV/FUSION neoantigen categories once the full Gide matrix exists.

2. **Riaz & Hugo RNA ↔ per-patient clonality** — the RNA-by-clonality
   interaction test set.
   - `data/registry/riaz_clonality.csv` — 68 patients, purity-free CCF proxy
     from the S3 MAF VAF; joins to the Riaz RNA cohort on the shared `Pt##`
     namespace. **51 patients carry RNA-seq + labeled response (11 R / 40 N).**
   - `data/registry/hugo_clonality.csv` — 38 patients, purity-corrected CCF
     (VAF·2/purity) from the S1D MAF + S1B purity; **21 R / 17 N**, 10 low-purity
     patients flagged. All 28 Hugo RNA runs join to it on base `Pt##`
     (Pt27A/Pt27B → Pt27).
   - **Caveat the modeling session must respect:** the two clonality tables are
     NOT cross-cohort comparable on *absolute* clonality (Riaz proxy median
     subclonal fraction 0.70 vs Hugo purity-corrected 0.25 — method, not
     biology). Use each within-cohort; for cross-cohort clonality, recompute both
     with one pipeline (PyClone-VI + FACETS/ASCAT purity).

## Known data limits

- **Gide has no WES**: clonality is inherently unavailable for that cohort, so
  the RNA-by-clonality test is Riaz + Hugo only.
- **No published FACETS purity in the Riaz supplements**: the Riaz clonality is
  a proxy by construction. A purity-corrected version needs raw exomes
  (open at ENA SRP095809, 177 WXS) run through a CCF caller, which is out of the
  current calls-only variant scope.

## Alignment output-format policy (durable — adopt for ALL alignment)

**Alignment output is COMPLETE, LOSSLESS, reference-based CRAM. Never a filtered BAM.**

Rationale: the project hypothesis names TE activation, fusion transcripts, and viral/repeat signal as
observable RNA phenotypes. Those live in unmapped and multi-mapped reads. An aligner run with
`--no-unal` / `-q 60` (as an AEI-only intermediate would use) discards exactly that signal — such a BAM
is a purpose-built intermediate, NOT a raw-read archive, and must not be reused for the raw-read / TE /
fusion branch.

Policy, implemented in `align_hisat2.sh`:
- Keep every read (no `--no-unal`, no MAPQ filter at align time).
- Emit coordinate-sorted, `.crai`-indexed CRAM with `--reference <genome.fa>`, default lossless quality
  (NO quality binning). Output dir `results/editing_crams/`.
- A complete CRAM is a full FASTQ replacement: `samtools fastq <cram>` regenerates the reads. Size is
  ~56% of FASTQ.gz (measured on the genuine complete CRAM ERR2208909: FASTQ.gz 5.54 GB → complete CRAM
  3.10 GB, 125.3M reads incl. 42M secondary — ~44% smaller and already aligned). A CRAM built from the
  filtered unique-read BAM is smaller (~1.5 GB) but is NOT complete — do not cite that as the ratio.
- Downstream tools filter at READ time, not align time: AEI applies `samtools mpileup -q 60` inside
  `compute_aei_fast.py --min-mapq 60`, which reads CRAM natively given the reference. Verified on a
  genuinely complete CRAM: ERR2208909 re-aligned from FASTQ with no filters (125.3M reads; 2.43M
  unmapped, 42.1M secondary, 59.1M MAPQ<60 — i.e. it contains exactly what a filtered BAM drops), then
  read-time `mpileup -q 60` reproduces the AEI from the old filtered-BAM alignment bit-for-bit
  (0.226002%, A>G=1040, A_cov=460172). So filtering at read time equals filtering at align time for AEI —
  the format switch loses nothing for AEI while preserving the raw reads for everything else.
- Do NOT auto-delete FASTQs in the align loop (the earlier version did, and removed 14/16 subset FASTQs).
  ENA is the canonical raw archive; `editing_subset_manifest.csv` holds accessions + `fastq_ftp` URLs and
  `prefetch_fastqs.sh` re-fetches with integrity checks, so re-materialization is a one-command,
  on-demand operation — but the pipeline should not throw raw data away as a side effect.

Aligner note: the `editing` conda env's STAR 2.7.11b arm64 build is broken (reads 0 input from every
FASTQ, confirmed on synthetic input) — use HISAT2 2.2.2 (the aligner that produced the repo's existing
BAMs), reusing `results/rnaseq_pilot_hisat2/genome/index/hisat2/`.

## Pointers

- Cohort rationale + exclusions: `docs/DATA_INVENTORY.md`, `data/registry/README.md`.
- Registry spreadsheet: `data/registry/DATASET_REGISTRY.xlsx`.
- Encoder decision (frozen-first, PEFT-gated): `docs/ENCODER_REVIEW.md`.
- Prior analysis handoff (the two directions these joins unblock):
  `HANDOFF_rna_state_next.md`.
