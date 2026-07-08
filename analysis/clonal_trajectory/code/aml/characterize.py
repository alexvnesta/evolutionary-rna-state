import numpy as np, pandas as pd, scanpy as sc, warnings
warnings.filterwarnings("ignore")
res=pd.read_csv("/tmp/aml/aml_clone_programs_classified.csv")
progs=["interferon_MHC","HSC_progenitor","myeloid_diff","cell_cycle"]
print("=== class means (raw program scores) ===")
print(res.groupby("class_k2")[progs].mean().round(3).to_string())
print("\n=== patients per class ===")
for c,g in res.groupby("class_k2"):
    print(f"class {c}: {g.patient.nunique()} patients, {len(g)} clones -> {sorted(g.patient.unique())}")
# z-scored class means to see defining axis
Z=(res[progs]-res[progs].mean())/res[progs].std()
Z["class_k2"]=res["class_k2"]
print("\n=== class means (z-scored) ===")
print(Z.groupby("class_k2")[progs].mean().round(3).to_string())

# genotyping concordance: fraction of mutant-transcript+ cells per clone
adata=sc.read_h5ad("/tmp/aml/aml_cnv.h5ad")
o=adata.obs
o["has_mut"]=o["MutTranscripts"].notna() & (o["MutTranscripts"].astype(str).str.strip()!="") & (o["MutTranscripts"].astype(str)!="nan")
mc=o[o.clone!="NA"]
print("\n=== mutant-transcript detection per clone (genotyping orthogonal check) ===")
tab=mc.groupby("clone").apply(lambda d: pd.Series({
  "n":len(d),"frac_mut_detected":round(d.has_mut.mean(),3)}))
print(tab.to_string())
print("cells with any mut transcript (clone cells):",int(mc.has_mut.sum()),"/",len(mc))
