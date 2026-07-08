import pandas as pd, numpy as np
cl=pd.read_csv("/tmp/aml/aml_clusters_all.csv")   # per-cluster cnv_ratio (pooled ref)
summ=pd.read_csv("/tmp/aml/aml_infercnv_summary.csv")
# strict per-clone gate: clone must be >=30 cells (is_clone True) AND cnv_ratio>1.5
cl["passes_strict_gate"]=cl["is_clone"] & (cl["cnv_ratio"]>1.5)
strict=cl.groupby("patient")["passes_strict_gate"].sum().rename("n_clones_strict_1p5x")
summ=summ.merge(strict,on="patient",how="left")
summ["multiclone_strict_1p5x"]=summ["n_clones_strict_1p5x"]>=2
summ.to_csv("/tmp/aml/aml_infercnv_summary.csv",index=False)
cl.to_csv("/tmp/aml/aml_clusters_all.csv",index=False)
print(summ[["patient","malig_cnv_ratio","n_clones","n_clones_strict_1p5x","multiclone_strict_1p5x"]].to_string())
print("\nmulticlone (Leiden, validation-based):", int((summ.n_clones>=2).sum()))
print("multiclone (strict per-clone >1.5x gate):", int(summ.multiclone_strict_1p5x.sum()))
print("patients with >=1 aneuploid clone:", int((summ.n_clones_strict_1p5x>=1).sum()))
