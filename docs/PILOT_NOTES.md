# Salmon pilot — environment notes (modeling session)

Independent lightweight de-novo quant arm (complements the sibling STAR/nf-core
pipeline). Local Mac, arm64, 64 GB RAM, no GPU.

## Toolchain (works on arm64)
- `salmon 2.3.1` (Rust rewrite / piscem) installs cleanly via bioconda on
  osx-arm64. **CLI changed from C++ salmon**: `salmon index` needs `--gencode`
  for GENCODE headers; no duplicate `-k`. `quant` keeps `--gcBias --seqBias
  --validateMappings -g <tx2gene>`.
- Reference: GENCODE v44 transcripts (252,835 tx). Transcriptome-only index =
  655 MB, builds in ~30 s with -p 8. (Decoy-aware/genome index is the sibling's
  STAR job.)

## The real bottleneck: ENA download bandwidth (~1.4 MB/s)
- Full FASTQ (~5.4 GB/sample) ⇒ ~14 h for 12 samples serial. Too slow.
- **Pilot fix: stream-subsample.** `curl <url> | zcat | head -n N*4 | gzip`
  takes the first N read pairs (mates stay in sync — FASTQ order preserved).
  3M pairs ≈ 190 MB/mate ≈ ~4 min/mate. Sufficient for gene-level TPM to prove
  de-novo signal + concordance. NREADS=0 in the script does full+MD5 for
  production. SIGPIPE from `head` is expected (guard with set +o pipefail).
- Peak disk bounded by stream-align-delete (delete FASTQ after quant).

## For production / full cohorts
Full-depth + MD5 verification, ideally on a host with real bandwidth and a
decoy-aware genome index. Consider remote compute if a cluster is added.
