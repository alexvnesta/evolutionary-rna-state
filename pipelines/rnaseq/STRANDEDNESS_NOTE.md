# Strandedness on Apple Silicon — READ BEFORE THE PRODUCTION RUN

nf-core/rnaseq's automatic strandedness detection (`strandedness=auto`) runs
`fq subsample` + Salmon on a read subset. **`fq` has no osx-arm64 conda build**,
so `auto` cannot run natively on this Mac (it fails at FASTQ_SUBSAMPLE_FQ_SALMON).

`samplesheet_selection.csv` (the 40-sample production sheet) is left at `auto`
to record the correct intent. Before the production run you MUST resolve
strandedness one of these ways:

1. **Derive it empirically (recommended, no guessing).** Salmon reports the
   observed library type in `results/.../salmon/<sample>/lib_format_counts.json`
   (`expected_format`). Run the pilot (which pseudo-aligns with Salmon on
   explicit lib type `A`=auto-detect *within Salmon*, which does NOT need fq),
   read the inferred strand from `lib_format_counts.json`, and set that value
   per cohort. Gide and Riaz are single-protocol within each cohort, so one
   determination per cohort suffices.
2. **Run the strandedness check once on a Linux/Docker machine** and copy the
   per-sample values back into the samplesheet.
3. If the library prep is known from the cohort's methods/GEO record, set it
   directly (Gide 2019 PRJEB23709 / Riaz 2017 GSE91061).

`samplesheet_pilot.csv` uses `reverse` ONLY to let the arm64 end-to-end
validation run proceed — this is a validation placeholder, NOT a verified
protocol value, and must not be propagated to production without confirmation.
