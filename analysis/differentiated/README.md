# differentiated — RNA-state-specific neoantigen modules

Modules that derive neoantigens a WES/annotation pipeline **cannot** see, from
de-novo transcriptomic phenotypes. Each derives candidate peptides from a
per-sample phenotype matrix and scores them through the **shared** antigen core
(`analysis/antigen_core`), so every `*_neoantigen_burden` is defined
identically and comparably.

## `splicing_neoantigen.py` — splicing-derived neoantigen burden

**Reference tool:** SNAF (Splicing Neo Antigen Finder), Li et al. 2024, *Sci.
Transl. Med.* (doi:10.1126/scitranslmed.ade2886; github.com/frankligy/SNAF).

**Feature produced:** `splice_neoantigen_burden` (int, per `run_accession` +
`cohort`).

### Pipeline

```
per-sample splice-junction COUNTS (STAR SJ.out.tab)
   │  call_neojunctions()      SNAF tumor-specificity gate
   ▼                           (count - normal_mean) >= t_min(20) AND normal_mean < n_max(3)
tumor-specific NEOJUNCTIONS
   │  retrieve_flanking_seq() + translate_junction()
   ▼                           SNAF 3-frame read-through of donor→acceptor exon
junction-SPANNING 8-11mer peptides
   │  count_binders()          SHARED MHCflurry engine (rank <= 2.0)
   ▼
splice_neoantigen_burden  (unique binding peptides)
```

```python
from analysis.differentiated.splicing_neoantigen import (
    parse_star_sj, splice_neoantigen_burden, build_feature_table)

junctions = parse_star_sj("SAMPLE_SJ.out.tab")            # STAR junction counts
alleles   = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", ...]  # from hla_typing
burden    = splice_neoantigen_burden(junctions, alleles) # -> int

# whole-cohort table (keyed run_accession + cohort):
tbl = build_feature_table(sample_junctions, sample_hla)  # -> DataFrame
```

### Install status & faithful-substitution note (arm64 macOS dev box)

`pip install snaf` is **not installable here**: SNAF 0.7.0 (and all versions
0.5.0–0.7.0) hard-pin `tensorflow==2.3.0`, which has no osx-arm64 / py3.11 wheel
(`pip install --dry-run snaf` → `ResolutionImpossible`). SNAF is also
Linux-oriented and needs an AltAnalyze exon-junction DB + a GTEx junction
control `.h5ad` not present in the sandbox.

Per the task's documented fallback, this module **wraps the SNAF algorithm
faithfully** rather than importing the unbuildable package. The algorithmic
core is a direct port of SNAF's source (`snaf/snaf.py`, `snaf/gtex.py`,
v0.7.0):

| SNAF source | this module |
|---|---|
| `crude_tumor_specificity` (gtex.py:252) | `is_neojunction` / `call_neojunctions` — identical count gate |
| `NeoJunction.in_silico_translation` (snaf.py:1033) | `translate_junction` — 3-frame (phase 0/1/2) read-through |
| `get_peptides` (snaf.py:1186) | `_get_peptides` — verbatim port: translate first-to-stop, continue into second, emit only junction-**spanning** k-mers |
| `subexon_tran`/`query_from_dict_fa` (AltAnalyze DB) | `retrieve_flanking_seq` — **SUBSTITUTION:** reads flanking exon sequence from GRCh38 FASTA + GENCODE annotation instead |
| `run_MHCflurry`/`run_netMHCpan` (binding.py) | `antigen_core.mhc_binding.count_binders` — **SHARED** engine |

So the peptide set is derived by SNAF's exact translation logic; only the
sequence-lookup and binding backends are swapped for this project's shared,
reproducible infrastructure. On the Linux pilot host SNAF *is* installable — if
the pilot swaps in upstream SNAF, only the `SPLICE_TRANSLATE` process body
changes; the channel contract and feature definition are unchanged.

### Batch robustness

The binder step is batch-invariant (MHCflurry percentile rank is
allele-calibrated against a fixed random-peptide background). The neojunction
**gate is count-based and therefore library-size sensitive** — the pilot sets
the tumor count threshold per platform and/or supplies a shared normal-junction
reference (recommend CPM-normalized counts). Report burden per clinical
context, not pooled across cohorts.

### Nextflow

`splicing_neoantigen.nf` is the subworkflow stub (mirrors `hla_typing.nf`):
`SPLICE_CALL_NEOJUNCTIONS → SPLICE_TRANSLATE → SPLICE_SCORE_BURDEN →
MERGE_SPLICE_BURDEN`, each wrapping a `splice_neoantigen_cli.py` subcommand so
the algorithm is single-sourced with the unit-tested Python. `stub{}` blocks
let `nextflow -stub-run` exercise the wiring without genome/MHCflurry present.

### Validation

`test_splicing_neoantigen.py` — synthetic genome only (never fabricated cohort
values). A designed junction, translated by SNAF's exact 3-frame logic, must
yield the influenza-M1 epitope **GILGFVFTL** (textbook HLA-A\*02:01 strong
binder) spanning the donor/acceptor boundary → burden ≥ 1; a low-count junction
must be dropped by the gate. Both `+` and `−` strand layouts are exercised.

```
cd <repo> && PYTHONPATH=. python analysis/differentiated/test_splicing_neoantigen.py
```

`DEMO_splice_neoantigen_burden_synthetic.csv` (in `results/features/`) is a
3-sample **synthetic** demonstration (sample ids `SYN_*`) — not real cohort
burdens.
