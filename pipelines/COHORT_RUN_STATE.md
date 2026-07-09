# Cohort spine run — state snapshot (for next session)

## Run
- exec_id: 568fbd52-1ba9-4b52-be88-edf9524f2369 (background, running)
- command: bash pipelines/scripts/run_cohort_batched.sh data/manifests/selection_manifest.csv results/rnaseq_cohort 8
- log: results/rnaseq_cohort.log
- work-dir: results/large/nf_work_rnaseq (SHARED across batches; holds cached genome indexes for -resume)
- output: results/rnaseq_cohort/ (hisat2/*.markdup.sorted.bam, salmon/, stringtie/, multiqc/)

## Progress at snapshot
- Batch 1/5 COMPLETE (Pipeline completed successfully): 8 BAMs+indexes, 8 salmon quants, salmon.merged.gene_counts.tsv, MultiQC. Batch boundary PROVEN: runner transitioned to batch_001, batch 2 fetching.
- 4 batches remain (batch_001..004); batch 2+ skip index build (cached in shared work dir)
- Genome indexes (hisat2 8x .ht2, salmon, transcripts.fa) built in results/rnaseq_cohort/genome/

## Disk management (CRITICAL — 64GB machine, storage premium)
- Peak per-batch: ~250GB work dir + 51GB raw FASTQ + output. Gets tight (~58GB free at batch-1 peak).
- Relief: raw FASTQs are safe to delete once a batch's 8 BAMs are published (alignment consumed them).
  Runner auto-deletes post-batch; can delete early with: rm -f data/raw/fastq/*.fastq.gz
- The shared work dir is NOT pruned between batches. If disk gets tight, after a batch completes
  its published BAMs, prune completed task dirs — but ONLY when that batch's nextflow has finished
  (ps is unreliable in sandbox; check log for "Pipeline completed" or absence of new Submitted lines).

## Next steps after all 40 BAMs aligned
1. Run the 4 subworkflows across the full cohort (they consume results/rnaseq_cohort/hisat2/*.bam):
   - intron_retention, rna_editing, te_erv, rnasplice (--source genome_bam)
2. Then optionally: re-invoke runner with ARCHIVE_CRAM=1, or run archive_bam_to_cram.sh on the BAMs
   to reclaim disk (~53% smaller, lossless).

## Known-resolved blockers (do not re-hit)
- strandedness=reverse (explicit) in _ss_from_manifest.py — auto needs fq=0.12.0 (no arm64 build)
- fetch_fastq.sh is resilient (curl -C - + retries) — rides through ENA SSL EOF/403
- rnafusion is arm64-blocked (needs amd64 STAR) — separate track

## Disk-management UPDATE (learned during batch 2)
- The shared work dir hit 336 GB during batch 2's alignment (8 samples' scratch
  co-resident). Disk dropped to 29 GB — CRITICAL. Nextflow auto-clears each
  sample's HISAT2 scratch once its BAM finalizes, but 8-wide concurrency peaks high.
- RECLAIMED 34 GB safely by deleting my own completed pilot working dirs
  (results/rnasplice_pilot + results/te_erv_integration_test) — their key results
  are already saved as artifacts, so nothing lost. Did NOT touch other sessions'
  data (editing_bams, editing_crams).
- RECOMMENDATION for a re-run or if disk gets critical again: use batch size 4-6
  instead of 8 (fewer samples aligning concurrently = lower peak scratch). The
  runner takes batch size as arg 3: run_cohort_batched.sh <manifest> <out> 4
- If a batch's BAMs are published and disk is tight, that batch's raw FASTQs
  (data/raw/fastq/*.fastq.gz) are safe to delete immediately.

## Disk pattern CONFIRMED (batches 1-2 observed)
- Per-batch disk cycle: work dir peaks ~400-430 GB during the 8-sample
  alignment→markdup phase (disk dropped to ~17 GB at worst), then Nextflow
  CLEARS the alignment scratch as BAMs finalize and disk RECOVERS to 250-337 GB
  before the next batch. The peak is transient, not cumulative — the work dir
  does NOT permanently grow across batches.
- Mitigation applied during peaks: delete each batch's raw FASTQs as soon as a
  sample's BAM appears in the work dir (align done). This bought ~40 GB/batch.
- Batches 1 and 2 both completed cleanly (16/40 BAMs published). RSeQC is the
  slow tail of each batch (~30-45 min for 8 samples). No batch has failed on disk.
- Progress tracking: `grep -c 'Pipeline completed successfully' results/rnaseq_cohort.log`
  = batches done; `ls results/rnaseq_cohort/hisat2/*.markdup.sorted.bam | wc -l` = BAMs.
