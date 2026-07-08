# analysis — analyses and benchmarks

Top-level analyses that sit on the `src/` modeling core, plus the baseline,
antigen-core, differentiated (per-phenotype), and pilot sub-packages.

## Modules
- `01_covariation.py` — covariation analysis testing falsifiable claim #1 (coordinated RNA-phenotype abnormalities).
- `baseline_benchmarks.py` — end-to-end baseline benchmark driver (TMB / SNV-indel neoantigen vs. RNA-state).

## Sub-packages
- `baseline/` — genomic baselines (standardized TMB, SNV/indel neoantigen).
- `antigen_core/` — shared antigen machinery (HLA typing, MHC binding).
- `differentiated/` — the five RNA-phenotype antigen arms.
- `pilot/` — de-novo Salmon pilot and its deepened follow-up.
