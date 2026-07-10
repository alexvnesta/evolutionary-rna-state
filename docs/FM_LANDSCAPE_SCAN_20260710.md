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
- Real long-context genomic FMs surfaced: Evo 2 (StripedHyena, primary), Caduceus (2403.03234, bi-directional equivariant), HyenaDNA (2306.15794, local dev stand-in), Nucleotide Transformer, and 2025 SSM-genomics work (2504.06304).
- **No genomic FM in this scan clearly supersedes Evo 2 as the long-context, single-nucleotide, autoregressive engine** needed for delta-likelihood scoring. Alternatives are mostly smaller / shorter-context.
- CAVEAT: arXiv relevance ranking is noisy (off-topic hits filtered out); this is a ~15-min scan, NOT a systematic benchmarked review. It supports 'Evo 2 is a defensible primary choice', not 'Evo 2 is provably SOTA'.
- Open question the protocol §6b resolves empirically: genome models are OOD on spliced/edited RNA, so an mRNA-specialised model may win for specific RNA-visible event classes. The ClinVar/splice-variant control task decides per class.