# Alternative-Splicing Tool Landscape for Bulk RNA-seq — 2026 Recommendation

**Scope:** 40 pre-treatment melanoma tumor RNA-seq samples (21 ICB responders / 19 non-responders), paired-end, reverse-stranded, ~40–80M reads/sample. Goal: event-level differential splicing (dPSI) as interpretable phenotype features for an "evolutionary RNA state" / tumor-antigenicity hypothesis, where **cryptic / novel junctions matter**. Compute: single Apple-Silicon (arm64) Mac, conda / no-Docker. Assets already in hand: Salmon transcript quantifications (all 40) + HISAT2 genome BAMs.

---

## 1. The current field (mid-2026)

Two 2025 reviews frame the state of play. The WIREs RNA overview (Tran et al., 2025) and the F1000Research practical-considerations review (Draper et al., 2025) both conclude there is **no single winner** — benchmark rankings flip depending on event type, sample size, annotation quality, and whether novel junctions are in scope. Draper et al. explicitly note that benchmarking studies "show no consensus on tool performance," and default to **DEXSeq and rMATS** on the strength of citation weight and active maintenance, while flagging **LeafCutter** as the go-to when annotation-free / novel-junction detection is needed.

The most quantitative recent head-to-head is Jiang et al. (2023, *Briefings in Bioinformatics*, event-level benchmark of 21 tools). Its findings drive most of the reasoning below:

- **At the event level, SUPPA ranked first, followed by LeafCutter, DARTS, and rMATS** — the four event-based tools outperformed gene-level approaches.
- **Novel-junction detection is a real weakness of annotation-anchored tools:** rMATS and DARTS *lost the ability to identify retained introns with novel junctions*; MAJIQ was the best at novel junctions but paid for it (true-positive rate dropped ~0.14, FDR rose ~0.15 at the gene level).
- **Tool unions beat single tools:** SUPPA+rMATS and SUPPA+LeafCutter were among the best-performing pairs; taking a union maximizes recall, taking an intersection tightens FDR.

FDR behaviour is the main knock on SUPPA2. In the 2023 event-level benchmark, annotation-anchored quant tools traded recall for precision differently than junction-based tools such as MAJIQ, which handled novel junctions best but at a true-positive-rate/FDR cost. SUPPA2's own paper (Trincado et al., 2018) reported FPR <5% on simulated data with ROC/PR area comparable to the alternatives it was tested against, so much of the field disagreement reflects defaults and thresholds rather than a fixed ranking. The consistent, defensible statement: **SUPPA2 is fast and accurate at moderate depth, but its detection counts and FDR drift more at larger sample sizes**, and it is annotation-bound. (Note: I did not retrieve a standalone MAJIQ-v2 re-benchmark with specific SUPPA2 FDR/FNR figures; the FDR characterization here rests on the Jiang 2023 benchmark and the SUPPA2 paper only.)

**Newer tool that has gained traction (2024–2026): SpliceWiz** (Wong et al., 2024, *Briefings in Bioinformatics*). It is a multi-threaded R/Bioconductor package that quantifies AS from BAM junction reads *and* intron-retention coverage, **leverages novel junction reads to detect cryptic splice sites and exons**, and runs GLM-based differential testing (DESeq2 / limma / DoubleExpSeq backends) optimized for large datasets. It directly targets the two things SUPPA2 cannot do (novel junctions + IR) while staying installable on a workstation.

---

## 2. Recommendation table

| Tool | Method basis | Fit for this design | arm64 / conda? | Verdict |
|---|---|---|---|---|
| **SUPPA2** | Event-level ΔPSI from transcript quant (annotation-derived events) | **Strong** — consumes existing Salmon quants directly (no realignment), fast, replicate-aware, top event-level ranking, interpretable per-event features. Blind to novel/cryptic junctions. | ✅ Pure-Python (numpy/pandas/scipy); native osx-arm64 | **PRIMARY (qualified)** |
| **SpliceWiz** | BAM junction + IR coverage, GLM (DESeq2/limma/DoubleExpSeq) | **Strong complement** — runs on the HISAT2 BAMs already in hand, detects cryptic splice sites/novel exons **and** intron retention, scales to large cohorts, interactive QC. Fills SUPPA2's novel-junction/IR blind spot locally. | ✅ Bioconductor R (compiles ompBAM C++ backend on osx-arm64) | **COMPLEMENTARY (primary cross-check)** |
| **LeafCutter** | Annotation-free intron-cluster excision, Dirichlet-multinomial GLM | Good for annotation-free **novel-junction confirmation**; used in cancer neoantigen-splicing pipelines. No native IR/event-type labels; R + regtools install is finicky. | ⚠️ regtools (C++) builds on arm64; R/differential side workable but fiddly | Optional 2nd cross-check |
| **rMATS-turbo** | Event-level counts, `--novelSS` for novel splice sites; scales to 10⁴ samples | Community-standard event tool; ideal on paper (novel-SS + scalable) and pairs well with SUPPA (top union). | ❌ **Blocked on osx-arm64** — `r-pairadise` conda-unsatisfiable; STAR realignment path also broke on arm64 | Run only on Linux/amd64 (remote) |
| **MAJIQ / MAJIQ v2** | Local splicing variations (LSVs), lowest reported FDR | Best FDR + novel-junction handling, but heavyweight; academic license + registration. | ❌ Linux-oriented C build; no supported osx-arm64 path | Linux/amd64 (remote) only |
| **DEXSeq + edgeR diffSplice** | Exon-bin differential usage (parametric) | Robust, well-maintained; your current fallback. Exon-level, not event-level ΔPSI; weaker for cryptic-junction interpretation. | ✅ Bioconductor R; native arm64 | Keep as sanity baseline |
| **Whippet** | Julia quasi-mapping over splice graph, node PSI | Fast, annotation+novel via contiguous splice graph. Lightly maintained; smaller community. | ⚠️ Julia runs on arm64; low recent activity | Niche / not recommended as primary |

---

## 3. Is SUPPA2 still defensible as primary in 2026? — **Qualified yes**

SUPPA2 has **not** been superseded for this design. It remains a top-ranked event-level tool (first at the event level in the 2023 21-tool benchmark), it is the only candidate that turns the **Salmon quants you already have** into event-level ΔPSI with essentially no extra compute, it is replicate-aware, and it is the one primary option that installs and runs natively on osx-arm64. For producing interpretable, named per-event dPSI phenotype features across 40 samples on a Mac, it is the right primary.

The qualification is important and maps exactly onto your hypothesis: **SUPPA2 derives its events from the reference transcriptome, so it is structurally blind to cryptic and novel junctions** — the very signal most relevant to tumor-specific antigenicity and neoantigen-generating aberrant splicing. It also shows more FDR/detection drift at larger sample sizes than junction-based tools. So SUPPA2 alone is insufficient for the cryptic-junction arm of the project — hence the mandatory complementary tool.

**Primary:** SUPPA2 (quant-based event ΔPSI). Use a stringent |ΔPSI| threshold (≥0.1) with the replicate-aware p-value, and treat its event catalog as annotation-anchored.

**Complementary:** SpliceWiz on the HISAT2 BAMs, for novel/cryptic splice sites + intron retention that SUPPA2 cannot see, and as an orthogonal-method concordance check on shared events. If you want a second, annotation-free novel-junction confirmation for the neoantigen arm, add LeafCutter (it is already established in cancer splicing-neoantigen pipelines). When a Linux/amd64 host becomes available, rMATS-turbo `--novelSS` is the natural third leg — the SUPPA+rMATS union was one of the best-performing combinations in benchmark, and rMATS-turbo scales cleanly to the full cohort.

---

## 4. arm64 / compute notes

- **Runs locally on the arm64 Mac (conda, no Docker):** SUPPA2 (pure Python), SpliceWiz + DEXSeq/edgeR (Bioconductor), LeafCutter (regtools compiles; R side workable). Junctions for SpliceWiz/LeafCutter come from the HISAT2 BAMs you already have — no STAR realignment needed.
- **Requires remote Linux/amd64:** rMATS-turbo (the `r-pairadise` dependency is conda-unsatisfiable on osx-arm64, and STAR realignment also fails on arm64) and MAJIQ (Linux C build, no supported osx-arm64 path; academic license). Flag both for a Linux box or cloud (rMATS-cloud CWL/WDL/Nextflow workflows exist) if you want to extend beyond the local combo.

---

## Bottom line

Keep **SUPPA2 as the primary** event-level ΔPSI caller — it is still a top-ranked event tool in 2026, uniquely turns your existing Salmon quants into interpretable per-event features with near-zero extra compute, and is the only strong primary that runs natively on the arm64 Mac. Pair it with **SpliceWiz** as the complementary cross-check: it runs on the HISAT2 BAMs you already have, adds the cryptic/novel-junction and intron-retention detection that SUPPA2 is structurally blind to (and that is central to the tumor-antigenicity hypothesis), and installs on arm64 via Bioconductor. SUPPA2 is defensible but not sufficient alone — the novel-junction gap is exactly your neoantigen signal, so the SUPPA2 + SpliceWiz pairing (optionally LeafCutter for annotation-free confirmation, and rMATS-turbo `--novelSS` when a Linux host is available) is the right local-first design.

---
*Key sources:* Draper et al. 2025 F1000Research (10.12688/f1000research.155223.2); Tran et al. 2025 WIREs RNA (10.1002/wrna.70030); Jiang et al. 2023 Brief Bioinform event-level benchmark (10.1093/bib/bbad121); Wong et al. 2024 SpliceWiz, Brief Bioinform (10.1093/bib/bbad468); Trincado et al. 2018 SUPPA2, Genome Biol (10.1186/s13059-018-1417-1); Wang et al. 2024 rMATS-turbo, Nat Protoc (10.1038/s41596-023-00944-2); Li et al. 2018 LeafCutter, Nat Genet (10.1038/s41588-017-0004-9).
