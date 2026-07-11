# Quantifying A-to-I RNA Editing from Bulk RNA-seq: 2026 Best-Practices Review and Recommendation

**Context:** 40 melanoma tumor RNA-seq (21 ICB responders / 19 non-responders), paired-end reverse-stranded, HISAT2 genome BAMs, MAPQ≥60 unique mappers. Goal: A-to-I editing burden as an interpretable phenotype feature (global index + regional/per-site calls; ADAR expression relevant). Platform: single Apple-Silicon (arm64) Mac, conda/no-Docker.

---

## 1. The two quantities and the established tools (2026)

A-to-I editing is read from RNA-seq as **A>G mismatches on the sense strand** (T>C on the antisense), because inosine is read as guanosine. Two complementary readouts:

- **Global Alu Editing Index (AEI)** — a single, robust, coverage-weighted ratio of edited (A>G) over total adenosine coverage pooled across *all* Alu adenosines genome-wide. It does not call individual sites, is largely insensitive to coverage/expression differences, and is the field-standard summary of ADAR activity per sample.
- **Per-site calling** — position-level editing sites and levels, enabling Alu-region editing maps, novel-site discovery, and hyper-editing.

**Landscape of tools:**

| Tool | What it does | Status / maintenance |
|---|---|---|
| **RNAEditingIndexer** (Roth, Ben-Aroya, Levanon; *Nat Methods* 2019, doi:10.1038/s41592-019-0610-9) | Reference implementation of the genome-wide AEI. C-based mpileup counter over a repeat (Alu) region set, with strand handling, SNP masking, read-end trimming built in. | Definitive AEI standard; repo stable but not actively developed. Docker-first, heavy native deps. |
| **REDItools2 / REDItools3** (Picardi/Pesole/Lo Giudice) | De-novo and known per-site A-to-I calling from BAMs; pileup-based with rich filters (base/mapping quality, homopolymer, read-position). REDItools3 is the current pure-Python/pysam generation, parallelized. | **Actively maintained**; REDItools3 on PyPI. Pairs with **REDIportal** (the A-to-I site database; latest update *NAR* 2023, doi:10.1093/nar/gkac1156). |
| **JACUSA2** (Piechotta/Dieterich; *Genome Biol* 2022, doi:10.1186/s13059-022-02676-0) | Java pileup engine; RNA-vs-RNA and condition contrasts, complex read signatures, replicate-aware; `JACUSA2helper` R package for filtering. Fastest tool in the 2023 benchmark. | **Actively maintained.** |
| **SPRINT** (Zhang; *Bioinformatics* 2017, doi:10.1093/bioinformatics/btx473) | SNP-free per-site + hyper-editing via SNV-duplet clustering. Highest REDIportal support and best SNP avoidance in benchmark, but fewest sites (cluster-constrained). | Maintained but dated; **its own aligner (BWA) and its own MAPQ-rewrite pre-step required — authors explicitly do not recommend it on splice-aware aligners; hyper-editing needs raw FASTQ.** |
| **RED-ML** (Xiong; *GigaScience* 2017) | Logistic-regression site classifier. | Human-only, tied to **GRCh37/BWA-junction reference**; limited. Poor fit. |
| **RES-Scanner2 / hyper-editing pipeline** (Porath/Levanon) | Hyper-editing detection by realigning unmapped reads after A→G transformation. | Niche; needs unmapped reads (see §5). |
| **Newer (2024–2026):** REDInet (*Brief Bioinform* 2025, doi:10.1093/bib/bbaf107) — TCN deep-learning site classifier trained on REDIportal; LoDEI (*2024*) — window-based **differential** editing; a 2025 cytoplasmic AEI-style index (bioRxiv). | ML/differential refinements layered on the pileup tools above. | Emerging; not yet standard. |

**Benchmark takeaway** (Morales, Rennie, Uchida, *BioTech* 2023, doi:10.3390/biotech12030056): across BWA/HISAT2/STAR, **REDItools2 recovered the most sites with strong REDIportal support** (recommended primary caller when runtime allows); **JACUSA2 was fastest and replicate-aware**; **SPRINT gave the most database-supported, SNP-clean calls but the fewest**. ~90% of REDIportal sites lie in Alu — which is exactly why an Alu-restricted index captures the dominant signal.

## 2. Is the samtools-mpileup-over-Alu-BED strategy legitimate and AEI-equivalent?

**Qualified yes — it *is* the AEI algorithm.** RNAEditingIndexer internally is a pileup counter over the Alu region set that pools A>G and A coverage into one weighted ratio. A single `samtools mpileup -l alu.bed` pass that tallies pooled A>G / (A+G) over Alu adenosines reproduces the same estimator and is ~10× faster than a per-interval pysam loop while staying arm64-native. It is methodologically sound provided the following are handled — each is a real pitfall:

- **Strand (critical for a reverse-stranded library).** A-to-I shows as A>G on the transcript's sense strand. Without strand assignment you conflate genuine sense-strand A>G with antisense T>C and mis-pool the ratio. Split/assign reads by inferred transcript strand (use the library's reverse-stranded flag + Alu strand) and count A>G on sense, T>C on antisense — do not naively sum both at every locus.
- **SNP masking.** Remove germline A/G (and T/C) polymorphisms (dbSNP common / gnomAD) before pooling, or true SNPs inflate the "edited" numerator. This is the single largest false-positive source without matched DNA.
- **Base/mapping quality.** Apply min base quality (≈Q25–30, stricter than mpileup's default Q13) and your MAPQ≥60 unique filter; **disable BAQ (`-B`)** — BAQ recomputation near mismatches systematically deflates editing counts.
- **Overlapping mate double-counting.** Paired-end mates overlapping at a site must be counted once. `samtools mpileup` performs overlap detection/quality-adjustment by default; confirm it is on and do not pass `--count-orphans`/`-A`.
- **Read-end and simple-repeat artifacts.** Trim the first/last few bases of reads (alignment-error prone) and restrict to Alu (not all SINEs/simple repeats), mirroring RNAEditingIndexer defaults.

Handled this way the strategy is a legitimate, publishable AEI. The one caveat: it is an *in-house re-implementation*, so it should be **validated against RNAEditingIndexer on a handful of samples** (run the reference tool via its Docker image on any Linux/amd64 host or a cloud VM — a one-time concordance check, not per-sample) to demonstrate numerical agreement before it carries the cohort.

## 3. Recommendation for our platform + goal

**Primary — optimized in-house single-pass samtools-mpileup AEI** (Alu BED, pooled A>G/(A+G), strand-aware, SNP-masked, BAQ off). It is arm64-native, fast, transparent, directly gives the interpretable global burden feature the hypothesis needs, and is defensible once validated against RNAEditingIndexer on a subset. The 200k-Alu subset is a reasonable speed fallback if needed, but the full-set single-pass is cheap enough (~minutes–low tens of minutes/sample) that subsetting is unnecessary and slightly weakens comparability to the published index.

**Complement — REDItools3** for per-site / Alu-region editing calls. It is **pip-installable (PyPI), pure-Python/pysam, arm64-native, actively maintained**, works directly on your HISAT2 BAMs, and connects to REDIportal for annotation. This supplies the regional/site-level layer (and novel-site discovery) that a global index cannot. Add **ADAR1/ADAR2 (and ADAR1-p150) expression** from your quantification as an orthogonal covariate — editing burden should track ADAR1-p150.

Keep **JACUSA2** in reserve as a second per-site caller for the responder-vs-non-responder *differential* contrast (replicate/group-aware, fast, Java/arm64-fine); **SPRINT** only if you later realign FASTQ for hyper-editing (it needs BWA + raw reads and is not recommended on HISAT2 BAMs).

## 4. arm64 / conda / Docker installability

| Tool | Install path | arm64 Mac / no-Docker? |
|---|---|---|
| In-house mpileup AEI | samtools + python (conda) | **Native — ideal.** |
| **REDItools3** | `pip install REDItools3` (pysam) | **Native — recommended complement.** |
| **JACUSA2** | `.jar` via conda-forge openjdk | **Native (JVM).** |
| SPRINT | conda/pip + BWA; dated deps | Installable but BWA-tied, FASTQ-based; awkward. |
| RED-ML | source; GRCh37/BWA-junction ref | Poor fit; skip. |
| **RNAEditingIndexer** | **Docker image (bedtools/samtools/Java native deps)** | **Docker-first, arm64-hostile — use only on a Linux/amd64 host or cloud for the one-time validation run.** |
| REDInet / LoDEI | pip / source | Native-ish; optional, not standard yet. |

## 5. Correctness caveats specific to AEI on HISAT2 (vs STAR) BAMs

- **MAPQ threshold.** HISAT2 unique mappers are **MAPQ 60**, not STAR's 255. Any tool/preset assuming 255 for "unique" would discard all your reads — your MAPQ≥60 filter is correct; just ensure the AEI/REDItools MAPQ parameter is set to 60, not a copied STAR default.
- **Aligner effect on counts.** In the 2023 benchmark HISAT2 yielded the **fewest** sites of the three aligners but with **good REDIportal support** — i.e., high precision, somewhat lower recall than STAR. For a *global index* this is fine (the ratio is stable); for absolute *site counts* be aware they are aligner-dependent, so keep the aligner fixed across all 40 samples for comparability and don't compare your counts to STAR-based literature counts.
- **Strand inference (A>G vs T>C).** With a reverse-stranded library, assign editing by transcript strand as in §2; mislabeled strand is the most common AEI error on stranded data.
- **Hyper-editing is invisible to this pipeline.** Reads with many closely-spaced edits often fail to map or fall below MAPQ 60 and are excluded from both the AEI and per-site calls. If hyper-editing burden is of interest, run the dedicated Levanon-style pipeline on **unmapped/low-MAPQ FASTQ** (A→G transformed realignment) as a separate track — it cannot be recovered from the MAPQ≥60 BAMs.
- **No matched DNA.** All 40 are tumor-only RNA; rely on dbSNP/gnomAD masking (and, for somatic melanoma, expect C>T UV signatures — do not mistake these for editing; only A>G/T>C over Alu counts).

---

## Bottom line

For an arm64/no-Docker platform with HISAT2 MAPQ≥60 BAMs, adopt the **optimized single-pass `samtools mpileup` over an Alu BED as the primary global AEI** — it is the reference AEI algorithm in a faster, native form, defensible once you confirm numerical concordance with RNAEditingIndexer on ~3–5 samples run on any Linux/amd64 or cloud host. Layer **REDItools3** (pip, arm64-native, actively maintained) on the same BAMs for per-site and Alu-region editing calls, and correlate both with **ADAR1-p150 expression**. Get strand assignment (reverse-stranded), SNP masking, BAQ-off, MAPQ=60, and mate-overlap handling right, and treat hyper-editing as a separate FASTQ-based track if needed. RNAEditingIndexer remains the gold standard but serves here as a one-time validation reference, not the production tool.
