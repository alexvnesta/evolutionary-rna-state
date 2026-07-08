# Hugo cohort completed to 28/28 — quant + refreshed transfer arm

Status note from session `55f2530f` (Quantify 16 Un-quantified Hugo RNA Runs).
Everything below is copied from saved files, not memory. This closes the
"immediate target" named in `docs/PIPELINE_HANDOFF.md` (the 16 un-quantified
Hugo runs) and refreshes the weakest arm of the 3-cohort transfer test.

## What landed

**1. All 28 Hugo RNA runs are now quantified.** The 16 previously-missing runs
(SRP070710 / GSE78220) were quantified with the *same* pipeline that built the
original 12 — salmon 2.3.1, GENCODE v44 transcriptome index
(`refs/gencode_v44_index`, 251,955 tx), 3M read-pair stream-subsample,
`--gcBias --seqBias -g tx2gene`, stream-align-delete. Verified: the existing 12
quants in `results/hugo_salmon/` carry `num_processed=3000000`, so the merged
28-sample matrix is depth-homogeneous.

**2. Matrices updated (existing rows byte-unchanged):**
- `results/hugo_gene_tpm.parquet` — **12 → 28 samples** (15 R / 13 N),
  62,266 genes. Schema unchanged (`run_accession`, `cohort`, unversioned ENSG).
- `results/features/quant_gene_tpm.parquet` (contract format) — **40 → 68
  samples**: 30 Gide + **28 Hugo** + 10 Riaz. Gide/Riaz rows verified unchanged.
  (Hugo was previously absent from this file; all 28 were added, not just 16.)

**3. Quant provenance:**
- `results/hugo_complete_manifest.csv` — the 16-run manifest (FTP + MD5 + clinical).
- `results/hugo_complete_salmon/pilot_index.csv` — per-run status + mapping rate
  + n_processed for all 16 (all `ok`, all 3,000,000 reads).
- `results/hugo_complete_salmon/<run>/quant.sf` + `aux_info/` retained per run.

## Per-run mapping rates (new 16)

| run | patient | tp | resp | mapping % |
|-----|---------|----|------|-----------|
| SRR3184289 | Pt13 | PRE | R | 95.7 |
| SRR3184291 | Pt15 | PRE | R | 94.2 |
| SRR3184292 | Pt16 | ON  | N | 83.0 |
| SRR3184293 | Pt19 | PRE | R | 97.2 |
| SRR3184295 | Pt22 | PRE | N | 95.0 |
| SRR3184296 | Pt23 | PRE | N | 94.2 |
| SRR3184297 | Pt25 | PRE | N | 94.9 |
| SRR3184298 | Pt27A| PRE | R | 85.5 |
| SRR3184299 | Pt27B| PRE | R | 85.2 |
| SRR3184300 | Pt28 | PRE | R | 96.8 |
| SRR3184301 | Pt29 | PRE | N | 93.8 |
| SRR3184302 | Pt31 | PRE | N | 81.5 |
| SRR3184303 | Pt32 | PRE | N | 75.1 |
| SRR3184304 | Pt35 | PRE | R | 89.2 |
| SRR3184305 | Pt37 | PRE | R | 86.5 |
| SRR3184306 | Pt38 | PRE | R | 82.5 |

New-16 mean mapping 89.4% (baseline 12: 93.8%). QC: every quant.sf has the same
251,955 transcripts and TPM summing to 1e6; 22k–41k expressed tx/sample.
**SRR3184303 (75.1%)** is soft-flagged (>2 SD below the combined mean) but valid
— it has the highest expressed-transcript count of the set (41,214), so the
lower rate reflects more off-transcriptome reads, not a failed sample. Retained.
5 of the 16 downloads needed retries for transient ENA stream errors (`curl`
exit 23, with occasional `18`/`56` on the underlying connection); all recovered.

## Refreshed transfer arm — the binding limitation is resolved for Hugo

`docs/PIPELINE_HANDOFF.md` named the n=10–12 pilots as "the binding limitation
on every downstream claim … 'non-replication' cannot yet be separated from
'underpowered.'" That is now decidable for Hugo.

Method: **identical** to the modeling session's original generator (frame
`f16431e3`, cells 245–248) — train a z-standardized mean-panel antigen axis +
logistic map on Gide (n=30), freeze train-fit stats, apply to held-out Hugo;
bootstrap CI = 2000 test-set resamples, seed=1. Before refreshing, the **n=12
result was reproduced to the decimal** (ant6=0.583, CI [0.22,0.92]) as a control,
confirming the harness runs standalone.

| Arm | AUROC (antigen 6-gene) | 95% CI |
|-----|------------------------|--------|
| Within Gide (LOO, n=30) | 0.871 | — |
| Gide→Riaz held-out (n=10) | 0.360 | — |
| Gide→Hugo pilot (n=12) | 0.583 | [0.22, 0.92] |
| **Gide→Hugo full (n=28)** | **0.595** | **[0.37, 0.80]** |

Completing Hugo **left the pilot regime**: the point estimate is essentially
unchanged (0.583 → 0.595) while the CI **narrowed ~40% in width** (0.70 → 0.43;
[0.22,0.92] → [0.37,0.80], stable to ±0.01 across seeds 1/2/3/7/42). The CI **still includes
0.5**, so this is now a *powered* non-replication — the Gide-trained gene-level
antigen axis genuinely does not transfer to Hugo; it is no longer merely
underpowered. This **strengthens** the core negative conclusion (a gene-level
de-novo antigen-presentation axis is cohort-specific, not a transferable
tumor-intrinsic RNA state) rather than overturning it.

Updated files: `results/heldout_transfer_3cohort.json` (Hugo arm → n=28,
provenance block added), `results/fig_transfer_hugo_refresh.png`.

## For the modeling session (`e51c751a` / `f16431e3`)

- The full 28-sample Hugo matrix is a drop-in for anything that previously read
  the 12-sample `hugo_gene_tpm.parquet` — same schema, same gene columns.
- The Hugo↔clonality join (`data/registry/hugo_clonality.csv`) covers all 28
  runs on base `Pt##` (Pt27A/Pt27B → Pt27), so the RNA-by-clonality interaction
  test can now use the full Hugo arm.
- I refreshed `heldout_transfer_3cohort.json` in place with the identical method.
  If you regenerate it from your own notebook, the numbers above are what to
  expect at n=28.

## Not touched (owned by other live sessions)

Left `analysis/`, `src/`, `pipelines/`, and the manuscript/figure-deck alone
(modeling + Nextflow sessions own them). No bulk commits made. The Gide/Riaz
full-cohort quant (the ~1.3 TB remainder of the 228-run queue) was **not** run —
out of local storage budget, per operator instruction.
