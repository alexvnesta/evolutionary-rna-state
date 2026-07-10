#!/usr/bin/env python3
"""Assemble the non-reference feature matrix from results/nonref_run/out/<acc>/ caller outputs.
Robust to partial outputs: editing (AEI) optional, TE/IR/splice from fast callers.
CPM-normalizes TE family counts within-sample. Writes results/predictor/nonref_matrix_cohort.parquet.
"""
import pandas as pd, numpy as np, os, glob

ERS = "/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"

def parse_fc(path):
    d = {}
    if not os.path.exists(path): return d
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 7: d[p[0]] = float(p[-1])
    return d

def assemble(out_root=f"{ERS}/results/nonref_run/out"):
    rows = []
    for sdir in sorted(glob.glob(f"{out_root}/*/")):
        run = os.path.basename(sdir.rstrip("/")); row = {"run_accession": run}; layers = []
        p = f"{sdir}/aei.tsv"
        if os.path.exists(p):
            a = pd.read_csv(p, sep="\t")
            row.update(editing_AEI_percent=float(a["AEI_percent"].iloc[0]), editing_SN=float(a["signal_to_noise"].iloc[0]),
                       editing_AG=float(a["AG_mismatches"].iloc[0]), editing_Acov=float(a["A_coverage"].iloc[0])); layers.append("editing")
        p = f"{sdir}/ir_summary.tsv"
        if os.path.exists(p):
            s = pd.read_csv(p, sep="\t"); ne = max(float(s["n_introns_evaluated"].iloc[0]), 1)
            row.update(ir_mean=float(s["mean_IR"].iloc[0]), ir_median=float(s["median_IR"].iloc[0]),
                       **{"ir_frac_gt0.1": float(s["n_IR_gt_0.1"].iloc[0])/ne}, ir_n_eval=ne); layers.append("ir")
        p = f"{sdir}/junctions.bed"
        if os.path.exists(p): row["splice_n_junctions"] = sum(1 for _ in open(p)); layers.append("splice")
        te = parse_fc(f"{sdir}/te_family.counts")
        if te:
            tot = sum(te.values()) or 1.0
            for fam, c in te.items(): row[f"te_{fam}"] = c/tot*1e6
            layers.append("te")
        row["_layers"] = "|".join(layers)
        if layers: rows.append(row)
    mat = pd.DataFrame(rows).set_index("run_accession")
    return mat

if __name__ == "__main__":
    mat = assemble()
    layers = mat.pop("_layers")
    mat.to_parquet(f"{ERS}/results/predictor/nonref_matrix_cohort.parquet")
    print(f"assembled {len(mat)} samples x {mat.shape[1]} features")
    print("layer completeness:"); print(layers.value_counts().to_string())
