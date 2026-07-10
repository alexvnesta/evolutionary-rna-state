# Genomic/RNA foundation-model landscape scan
Date: 2026-07-10 | Source: arXiv (via literature MCP) | Purpose: evidence the Evo 2 model choice (protocol §6b)

## Query hits (genomics-filtered)

**Caduceus DNA language model**
- `2403.03234` — Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling
- `q-bio/0703003` — Effect of Internal Viscosity on Brownian Dynamics of DNA Molecules in Shear Flow
- `2105.03431` — DNA Nanotechnology Meets Nanophotonics
- `1607.00266` — The Art of DNA Strings: Sixteen Years of DNA Coding Theory

**Nucleotide Transformer genomics**
- `q-bio/0702036` — Nucleotide Distribution Patterns in Insect Genomes
- `2504.06304` — Leveraging State Space Models in Long Range Genomics
- `2505.08918` — When repeats drive the vocabulary: a Byte-Pair Encoding analysis of T2T primate genomes
- `1208.0133` — On pairwise distances and median score of three genomes under DCJ
- `2306.15794` — HyenaDNA: Long-Range Genomic Sequence Modeling at Single Nucleotide Resolution

**HyenaDNA long-range genomic**
- `2306.15794` — HyenaDNA: Long-Range Genomic Sequence Modeling at Single Nucleotide Resolution
- `2504.06304` — Leveraging State Space Models in Long Range Genomics
- `2201.08443` — Diversifying the Genomic Data Science Research Community

**Evo genome foundation model DNA**
- `2205.05897` — CAGI, the Critical Assessment of Genome Interpretation, establishes progress and prospects for computational genetic variant interpretation methods
- `q-bio/0703003` — Effect of Internal Viscosity on Brownian Dynamics of DNA Molecules in Shear Flow
- `2603.27465` — Poisoning the Genome: Targeted Backdoor Attacks on DNA Foundation Models
- `2603.06950` — How Private Are DNA Embeddings? Inverting Foundation Model Representations of Genomic Sequences
- `1607.00266` — The Art of DNA Strings: Sixteen Years of DNA Coding Theory

**RNA foundation model splicing isoform**
- `1502.05667` — Towards de novo RNA 3D structure prediction
- `2511.02622` — Machine Learning for RNA Secondary Structure Prediction: a review of current methods and challenges
- `1304.5952` — Methods to study splicing from high-throughput RNA Sequencing data
- `1602.06317` — Statistical modeling of isoform splicing dynamics from RNA-seq time series data
- `1603.05915` — MSIQ: Joint Modeling of Multiple RNA-seq Samples for Accurate Isoform Quantification

## Interpretation
- Genomic FMs the queries ACTUALLY returned (with the arXiv id located): Caduceus (2403.03234, bi-directional
  equivariant), HyenaDNA (2306.15794, local dev stand-in), and 2025 SSM-genomics work (2504.06304); Nucleotide
  Transformer appeared only as a keyword, not a located record.
- **Evo 2 / StripedHyena was NOT located by any query** (the "Evo genome foundation model DNA" query returned
  only off-topic hits — CAGI, DNA Brownian dynamics, an embedding-inversion-attack paper, a backdoor-poisoning
  paper, a DNA-coding-theory paper). So this scan does **not** independently confirm Evo 2's status one way or
  the other — it only shows the *alternatives* that surfaced are mostly smaller / shorter-context than Evo 2 is
  known (from its 2025 publication) to be. Evo 2's own capabilities here rest on prior knowledge, not this scan.
- What the scan supports: **no returned alternative is an obvious long-context superset of Evo 2**; it does NOT
  establish "nothing supersedes Evo 2", because a superior model absent from these noisy keyword queries would
  not appear. Treat as "no red flag found in a quick scan", not "confirmed non-supersession".
- CAVEAT: arXiv relevance ranking is noisy (off-topic hits filtered out); this is a ~15-min scan, NOT a
  systematic benchmarked review. It supports 'Evo 2 is a defensible primary choice', not 'Evo 2 is provably SOTA'.
- This is exactly why model choice is NOT settled by this scan: the protocol's §6b ClinVar/splice control
  benchmark is the actual gate that decides the model (per event class), empirically, before any GPU spend.