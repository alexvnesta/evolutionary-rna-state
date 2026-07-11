# TE/ERV RNA-seq quantification — 2026 best-practice review and recommendation

**Context:** 40 melanoma tumor RNA-seq samples (21 ICB responders / 19 non-responders), paired-end, reverse-stranded. On disk: HISAT2 genome BAMs filtered to **unique mappers** (MAPQ ≥ 60, `-F 256`) + Salmon transcript quants + GENCODE v46/GRCh38. **FASTQs deleted** (re-fetchable from SRA/AWS, ~127 GB + hours). Current `te_erv` subworkflow = Telescope + TEtranscripts, both requiring a `bowtie2 -k 100` multi-mapping realignment from FASTQ. Goal: TE/ERV (esp. ERV/LTR) expression as an interpretable phenotype feature for a tumor "evolutionary RNA state" / antigenicity hypothesis; locus-level ERV desirable, subfamily acceptable. Platform: single arm64 Mac, conda/no-Docker.

---

## 1. The central problem: young/active ERVs live in the multi-mappers

Every credible TE tool exists to solve one problem — reads from repetitive elements map to many loci. The direction of the bias is the crux for this project. The TEtranscripts authors state it plainly: using only uniquely-mapped TE reads **"will bias read assignment away from the youngest TE sub-families and toward the older... subfamilies with higher uniquely mappable content"** — because active/young elements have accumulated few polymorphisms and are therefore represented almost entirely by multi-mappers (Jin et al. 2015, Bioinformatics, doi:10.1093/bioinformatics/btv422). Independent work confirms unimappers are enriched for *old* TEs and multimappers for *young* ones (e.g. AluYa5, L1HS) at p < 2.2×10⁻¹⁶ (Berrens/Corbo-style analyses; biorxiv 2023.07.04.547702).

**Why this matters specifically for you:** the antigen-relevant ERVs — young HML-2/HERV-K proviruses, HERVH, young L1 — are exactly the ones that a unique-mapper-only count undercounts most. So a unique-only feature is not just "noisier," it is *systematically depleted in the biologically interesting direction*.

## 2. Can any tool quantify TE/ERV from the EXISTING unique-mapper BAMs without realignment?

Yes — but only at the **subfamily/family level, as a conservative (biased-low) proxy**. The multi-mapping EM tools cannot help you here because the multi-mapping alignments have already been discarded from your BAM (`-F 256` removed secondaries; MAPQ ≥ 60 kept only uniques). No tool can re-derive alignments that are not in the file, and `samtools fastq` on the filtered BAM would only recover the unique subset — not the multi-mappers you actually need. What can run on your BAMs as-is:

- **TEtranscripts `--mode uniq` / TEcount** — accepts any BAM, counts only unique reads → **subfamily-level** TE counts. On a BAM that is already unique-only, this is the natural mode. This is literally the "unique counts" comparator used as the low-resolution baseline in the Telescope paper (htseq-count over the same loci). Defensible as a floor; not as truth.
- **TElocal `--mode uniq`** — locus-level relative of TEtranscripts; also accepts a BAM. But it ships pre-built `.locInd` indices only for standard genome builds, its EM design assumes a multi-mapping BAM, and in `uniq` mode it inherits the same young-ERV depletion.
- **featureCounts over a RepeatMasker/HERV GTF** — simplest option; by default counts only uniquely-mapped reads, and your reverse-stranded protocol lets it separate genic from TE signal. A locus-level benchmark (Schwarz et al. 2022, Front Genet, doi:10.3389/fgene.2022.1026847) found featureCounts unique/fraction/random modes and even the best tool (TElocal) all produce a locus-level false-discovery count exceeding the number of correctly recovered active loci — i.e. **locus-level unique counting is not trustworthy**; subfamily-level aggregation is where it holds up.

**Verdict on unique-mapper-only quant:** *qualified.* Scientifically defensible at the **subfamily/family level** as a conservative lower-bound feature (and honestly, if reported as such, publishable). **Not** defensible at locus level, and **not** an unbiased readout of the young/active ERVs your antigenicity hypothesis centers on. What is lost: the youngest LTR/ERV and L1 families — precisely the antigen-relevant fraction — plus all locus identity.

## 3. Salmon-based routes (SalmonTE / REdiscoverTE) — cheaper, but still need FASTQ

- **SalmonTE** — Salmon quant against a Repbase TE-consensus transcriptome; ~20× faster than TEtranscripts. **Subfamily-level only**, age-related biases, and in the 2025 benchmark it was the *worst* at quantifying young families. **Needs FASTQ** (Salmon reads FASTQ, not your genome BAM).
- **REdiscoverTE** — Salmon against a whole-transcriptome index (GENCODE transcripts + introns + >5M RepeatMasker sequences); the standard for TCGA-scale TE-immune/antigenicity work (Kong et al. 2019, Nat Commun, doi:10.1038/s41467-019-13035-2). More comprehensive than SalmonTE, gives locus-resolved RepeatMasker output. **Needs FASTQ.**

Neither runs from your Salmon *gene/transcript* quants — they require a Salmon index built against a TE/repeat reference and a re-quant from reads. So a Salmon-based path is **cheaper than bowtie2 -k 100** (quasi-mapping vs. very-sensitive-local -k 100 alignment) but **still requires the FASTQ refetch**. Your existing Salmon quants against GENCODE transcripts do not contain the repeat signal.

## 4. The 2025/2026 landscape and where the field actually sits

The most authoritative current source is a **September 2025 benchmark of 16 teRNA tools** across 180 datasets (She, Wang & Yang, bioRxiv doi:10.1101/2025.09.30.679421). Headline results relevant to you:

- **Telescope was the most accurate and robust quantifier at both family- and unit(locus)-level.** It remains the reference standard for locus-level ERV work — and it is what a May 2026 acral-melanoma HERV paper and a 2025 HCC HERV paper both used (`bowtie2 --very-sensitive-local -k 100 --score-min L,0,1.6` → Telescope `retro.hg38.v1`, 1,054 named HERV loci). Your existing choice is well-founded.
- **RSEM / Kallisto win at exon- and transcript-level**; the benchmark argues the **TE-exon level** is the accuracy/resolution sweet spot, and packages this as a pipeline (**TERA**) — but its detect module needs ~16 h / 100 GB RAM / FASTQ, so it is heavier than what you have, not lighter.
- **Every accurate quantifier needs a multi-mapping alignment or FASTQ.** Telescope, TEtranscripts `multi`, TElocal `multi`, SQuIRE, RSEM, Kallisto, SalmonTE, REdiscoverTE — none accepts a unique-only BAM as a route to unbiased young-ERV signal.

**Newer 2026 tools worth knowing:**
- **MAJEC** (Calico; bioRxiv 2026, title "unified gene, isoform, and locus-level transposable element quantification from RNA-seq", biorxiv.org/content/10.64898/2026.04.10.717472 — GitHub calico/majec) — per its title, a unified method giving **gene, isoform, and locus-level TE quantification in one pass**. I located its listing/abstract but did not retrieve its full text, so I am *not* asserting its quantitative accuracy claims here; treat it as an emerging locus-level alternative to verify directly before adoption. Like the other accurate locus tools it operates on a multi-mapping alignment (not a unique-only BAM), so it would still require a realignment from FASTQ; if that realignment uses STAR it is faster than `bowtie2 --very-sensitive-local -k 100`.
- **LATTE** (Mar 2026) and **LocusMasterTE** (2025) — locus-level EM tools; LocusMasterTE needs paired long-read TPM, so not applicable.
- **ERVmap** — BWA + very stringent unique-best-match filtering over 3,220 proviral ERV loci. It is essentially a curated unique-mapper approach (over-discards ambiguous reads) and **still realigns from FASTQ with BWA**; its annotation is far smaller than Telescope's (~15k elements). Not a BAM-reuse shortcut.

## Recommendation table

| Tool | Needs multimap realign / FASTQ? | Works on your unique BAMs as-is? | Level | arm64 / conda | Verdict for this project |
|---|---|---|---|---|---|
| **TEtranscripts `--mode uniq` / TEcount** | No | **Yes** | Subfamily | Yes (bioconda, pysam) | **Interim floor feature** — biased low on young ERVs; label as lower bound |
| **featureCounts over RepeatMasker/HERV GTF** | No | **Yes** | Subfamily (locus unreliable) | Yes (subread, bioconda) | Simplest no-refetch option; use reverse-strand; subfamily only |
| **TElocal `--mode uniq`** | No (but designed for multi) | Yes (needs prebuilt locInd) | Locus (unreliable in uniq) | Yes (bioconda) | Locus counts untrustworthy from unique-only; limited value |
| **Telescope** | **Yes** (bowtie2/STAR -k 100 + FASTQ) | No | Locus (HERV) | Yes (bioconda, pysam/scipy) | **Definitive locus-ERV**; your validated choice; 2025 benchmark #1 |
| **MAJEC (2026)** | **Yes** (multi-mapping BAM + FASTQ) | No | Gene+isoform+locus TE | Likely (Python; verify) | Emerging locus alternative (full text not retrieved — verify claims); STAR realign faster than bowtie2 |
| **REdiscoverTE** | **Yes** (Salmon + FASTQ) | No | Locus/subfamily (RepeatMasker) | Salmon arm64 ok; R | Best for TCGA-style TE-antigenicity; still needs FASTQ |
| **SalmonTE** | **Yes** (Salmon + FASTQ) | No | Subfamily | Salmon arm64 ok; Py wrapper | Fast but worst on young families; skip |
| **SQuIRE** | **Yes** (STAR + FASTQ) | No | Locus | Linux/older Py2 deps — flag | Dated deps; not worth it on arm64 |
| **ERVmap** | **Yes** (BWA + FASTQ) | No | Locus (3,220 ERVs) | bwa arm64 ok; perl/py | Small annotation; over-discards; no BAM shortcut |

*arm64 note:* the Python/pysam tools (TEtranscripts, TElocal, Telescope, MAJEC), subread/featureCounts, bowtie2, bwa and Salmon all have osx-arm64 conda builds or build cleanly. **STAR** on osx-arm64 is the one to verify — bioconda coverage is inconsistent and it may need Rosetta/x86 emulation; if you adopt MAJEC or want STAR realignment, test the STAR install first or fall back to bowtie2 `-k 100` (which Telescope already uses and which is arm64-clean). ERV annotations are freely available: Telescope `retro.hg38.v1` (1,054 HERV loci), RepeatMasker/Dfam, ERVmap's 3,220-locus set, gEVE.

## Bottom line

You can avoid the refetch **only for a conservative, subfamily-level TE/ERV feature**: run TEtranscripts `--mode uniq` (or featureCounts over a RepeatMasker+HERV GTF) directly on your existing HISAT2 unique-mapper BAMs — zero download, minutes of compute, arm64-clean — and report it explicitly as a biased-low floor that under-represents the young/active ERVs. That is defensible and may be enough if the ERV feature is one of several inputs to a phenotype model. **It is not defensible as locus-level ERV quantification, and it structurally undercounts exactly the young HERV-K/HERVH/L1 elements your antigenicity hypothesis cares about.** For the locus-resolved ERV readout the hypothesis actually needs, there is no shortcut around a multi-mapping realignment, and since the FASTQs are gone that means the ~127 GB refetch — but you are not locked into bowtie2: keep Telescope (still the 2025 benchmark's best locus/family quantifier and your validated path), or realign once with STAR `-k 100` and feed Telescope (MAJEC, 2026, is an emerging locus-level alternative to evaluate directly before adoption). Recommended plan: ship the unique-BAM subfamily feature now as an interim/floor, and schedule the refetch + `-k 100` + Telescope run to produce the definitive locus-level ERV feature.
