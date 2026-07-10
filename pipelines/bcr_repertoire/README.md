# pipelines/bcr_repertoire — BCR/TCR repertoire reconstruction (TRUST4)

Raw-read arm that reconstructs immunoglobulin (and TCR) sequences from bulk
RNA-seq, feeding the somatic-hypermutation / class-switch / clonality features
in `analysis/differentiated/bcr_shm.py`.

## What it does

`run_trust4_pilot.sh` mirrors the salmon pilot's **stream-align-delete**
pattern: for each run in the manifest it streams read-pairs from the AWS SRA
mirror / ENA, runs TRUST4, keeps only the small report files
(`*_cdr3.out`, `*_report.tsv`, `*_airr.tsv`), and deletes the FASTQ. Peak disk
is ~1 sample.

## Tool

TRUST4 (Song et al., *Nat Methods* 2021, doi:10.1038/s41592-021-01142-2) is
built from source under `tools/TRUST4-1.1.5` (arm64-clean — the osx-arm64 conda
build is unavailable). References: `tools/hg38_bcrtcr.fa`,
`tools/human_IMGT+C.fa`. `tools/` is git-ignored (vendored binaries); rebuild
from the release tarball if absent.

## Downstream

```
run_trust4_pilot.sh  ->  <OUT>/<sample>_{cdr3.out,report.tsv,airr.tsv}
        │
        ▼  analysis/differentiated/bcr_shm.py :: build_bcr_features(OUT)
per-sample BCR features (SHM rate, isotype fractions, clonality, diversity)
```

Nextflow form of the same flow: `analysis/differentiated/bcr_shm.nf` with
`modules/{trust4_assemble,bcr_shm_features,merge_bcr_features}.nf`.

**Status:** pilot script; part of the differentiated (RNA-native) feature axis.
Not yet run at cohort scale. The canonical cohort non-reference feature build is
`results/nonref_run/` (editing / TE / IR / splicing); BCR/SHM is an additional
differentiated arm layered on the same aligned reads.
