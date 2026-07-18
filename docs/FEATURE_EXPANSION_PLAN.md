# Non-reference feature expansion — inventory, gaps, and build order

_Coordination session (`837512d2`). Grounded on the actual tested matrix
(`nonref_matrix_cohort.parquet`, 66 features) and repo caller inventory, not memory._

## 1. What is actually measured today (the tested block = 66 features)

| Layer | Features | Caller | In powered test? |
|---|---:|---|---|
| TE / ERV families | **57** | Telescope / TEtranscripts, family-level | ✅ |
| Intron retention | 4 | featureCounts intron/exon → IR ratio | ✅ |
| RNA editing (ADAR, A→I) | 4 | Alu Editing Index (`AEI_percent, SN, AG, Acov`) | ✅ |
| Splice junctions | 1 | raw novel-junction count | ✅ |

**86% of the non-ref matrix is family-level TE.** The immune floor it is tested
against is ~7 features (`gep_tcell_inflamed, ifng_score, teff, tgfb,
teff_tgfb_balance, tmb, hla_het`). Consequence for interpretation: the n=106
two-block null is a null on *aggregate TE/editing/IR burden*, NOT on "the
tumor's antigenic RNA state" — the antigenicity layers are barely instrumented.

## 2. Exists as code but NOT run at cohort scale (demo/synthetic only)

The four **antigen** callers below have a caller + a test, but the tests are
**synthetic-fixture logic checks** (peptide generation + MHC binding logic is
validated; no real per-sample burdens computed). **None of the four antigen
callers is in the matrix.** (The two burden callers `intron_retention.py` and
`rna_editing.py` DO already feed the matrix — see §1 table — and are listed here
only for completeness of the caller inventory.)

- `fusion_antigen.py` (.nf) — Arriba-based fusion detection + fusion-junction neoepitopes. Closest to ready. NOT in matrix.
- `splicing_neoantigen.py` (.nf) — SNAF-style neojunction peptides. NOT in matrix.
- `snv_indel_neoantigen.py` — SNV/indel neoantigens (missense + frameshift). NOT in matrix.
- `te_antigen.py` (.nf) — TE/ERV-derived peptides (needs locus-level TE). NOT in matrix.
- `intron_retention.py`, `rna_editing.py` — burden callers, cohort versions ALREADY feed the matrix (in powered test).

Dependency gap for all neoantigen arms: **HLA typing exists for only 16/106
samples** (`results/hla/`, the editing subset). Cohort neoantigen burden needs
HLA class-I typing on all 106 first (`hla_typing/` scaffolding present).

## 3. Not present at all (named in hypothesis or standard, no caller)

- **Cryptic / novel ORFs** — explicitly in the project hypothesis; no caller.
- **circRNA** (back-splice junctions) — distinct non-ref class; not the linear splice count.
- **APOBEC (C→U) editing** — we only measure A→I (ADAR). Immune-relevant, separate signal.
- **Locus-level TE** — all 57 TE features are family-level; locus Telescope designed but never verified on real data. Specific ERV loci (not families) are where antigen signal lives.
- **BCR/TCR repertoire** — TRUST4 exists but ran on 2 pilot samples only.
- **Junction annotation layer** (see §4).
- Lower priority: APA / 3′UTR shortening, allele-specific expression.

## 4. Junction support / long-read annotation (user question)

Today a splice junction is "novel" purely by **absence from GENCODE v46** — a
binary test. Upgrade it to a graded, evidence-stratified label:

- **GENCODE TSL (1–5)** of flanking annotated junctions: a novel junction is
  more confidently tumor-specific when neighbors are TSL-1 (well-supported),
  less so in TSL-4/5 regions where "novel" may just mean "under-annotated."
- **Long-read isoform support**: we have NO long-read data on these ICB
  patients (all short-read cohorts), so this is a **reference-catalog
  intersection**, not new sequencing — intersect novel junctions against a
  long-read SJ catalog (GTEx/ENCODE PacBio/ONT) and assign **SQANTI**
  categories (FSM/ISM/NIC/NNC). Reclassifies each junction as known /
  novel-in-catalog / genuinely unannotated using orthogonal evidence.

**Verdict:** this is an *annotation/filtering step applied at junction-calling
time* (makes the "novel" label specific), NOT a downstream step. It upgrades the
single `splice_n_junctions` count into confidence-stratified counts.

## 5. The unified Nextflow pipeline — what it does and doesn't wire

`pipelines/main.nf` (`unified_main.nf`) is a real DSL2 wrapper, but it wires
only **3 arms**, all consuming the shared spine BAM:

- ✅ RNA_EDITING, ✅ INTRON_RETENTION, ✅ TE_ERV (locus-level behind `--te_locus`).

Its own header states splicing (rnasplice) and fusion (rnafusion) are run as
**separate nf-core entries** (they re-align from FASTQ, not the spine BAM) and
merged later by the matrix builder (`build_powered_n106_matrix.py`; the pipeline
header refers to a `build_nonref_matrix.py`, which is not present in the repo —
the powered builder is the actual merge script). Neoantigen / ORF / circRNA / APOBEC /
repertoire / junction-annotation are **not in any pipeline** — standalone
callers or absent.

So "the Nextflow pipeline does it all" is true only for the 3 tumor-intrinsic
burden layers. The antigenicity half was never integrated — and that is the
half that connects the RNA layer to the ICB/antigen endpoint.

## 6. Build order (decided: ready-caller wave first, dev in parallel)

Substrate advantage: the box is producing a **uniform STAR re-alignment of all
106** (currently ~93/106). One aligner, one BAM per sample = the clean substrate
a unified feature pipeline needs (no aligner confound, no per-arm re-align).

**Wave A — run existing validated callers on the 106 STAR BAMs (once aligned):**
fusion (Arriba), ADAR editing, intron retention, TE family + **locus-level**,
splice-junction calling. Produces real cohort features from code that already
works. De-risks the powered test soonest.

**Wave B — develop the missing antigenicity arms (in parallel with A):**
1. HLA class-I typing on all 106 (blocks all neoantigen arms).
2. Wire the three neoantigen callers (SNV/indel, splice, TE) to real BAMs + HLA.
3. New callers: circRNA (CIRCexplorer/CIRI2), cryptic ORF, APOBEC editing.
4. TRUST4 BCR/TCR on all 106.
5. Junction-annotation layer (§4): GENCODE TSL + SQANTI vs long-read SJ catalog.

**Wave C — integrate every arm into `pipelines/main.nf` and merge** into the
powered n=106 matrix via `build_powered_n106_matrix.py`. Re-run the two-block
test with the full antigenicity block, not just TE burden.

All heavy compute runs on the Ubuntu box (24 cores + GPU, 434 GB free after the
disk relief). Alignment ownership stays with the Mac session; feature callers
consume the published BAMs read-only until alignment completes.
