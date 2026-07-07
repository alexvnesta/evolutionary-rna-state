# Published supplemental data

Author-provided supplemental tables from the two melanoma ICB papers whose raw
RNA reads we use, obtained during the hackathon and **verified against the
published papers** (see per-cohort READMEs). These give per-sample WES,
neoantigen, clonality, and signature data at finer granularity than the
harmonized cBioPortal release.

## Layout

```
supplements/
  hugo2016/
    raw/     original publisher files (mmc1.xls, mmc2.xlsx, mmc3.pdf) + checksums.json
    clean/   extracted, cleaned CSVs (one per sheet) + manifests
  riaz2017/
    raw/     original publisher files (mmc1-6.xlsx) + checksums.json
    clean/   extracted, cleaned CSVs (one per sheet) + manifests
```

**Git policy:** the large raw mutation workbooks (Hugo `mmc1.xls` 15 MB) are
git-ignored; the cleaned CSVs and the small raw files are committed. All raw
files are preserved locally with MD5s recorded in each `raw/checksums.json`.

## Provenance & verification

- **Hugo et al. 2016**, *Cell* 165:35–44 (PMID 26997480, GEO GSE78220).
  Verified: WES median coverage 140× (`AvgCov` med 140.3), median 489
  non-synonymous mutations (`TotalNonSyn` med 489, n=38), BRCA2-responder
  column present, IPRES gene sets present. See `hugo2016/README.md`.
- **Riaz et al. 2017**, *Cell* 171:934–949 (PMID 29033130). Verified: 73
  clinical patients, mutation-load median 182 (paper 183; range 1–7360 exact),
  68 patients with somatic calls, longitudinal genomic-availability matrix.
  See `riaz2017/README.md`.

Files supplied by the user (public author-provided supplements). Public data
only; consistent with the compliance boundary (data, not code).
