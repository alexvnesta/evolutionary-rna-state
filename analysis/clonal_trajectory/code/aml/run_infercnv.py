import os, warnings, numpy as np, pandas as pd, scanpy as sc
np.random.seed(0)
warnings.filterwarnings("ignore")
import infercnvpy as cnv
# monkeypatch process_map -> serial (numba/parallel discipline)
import infercnvpy.tl._infercnv as _icnv
def _serial_map(fn, *iterables, **kw):
    return [fn(*args) for args in zip(*iterables)]
_icnv.process_map = _serial_map

adata = sc.read_h5ad("/tmp/aml/aml_raw.h5ad")
adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# gene positions
gp = pd.read_csv("/tmp/cscc/gene_pos.csv")
gp.columns=["symbol","chrom","start","end"]
gp = gp.drop_duplicates("symbol").set_index("symbol")
common = adata.var_names.intersection(gp.index)
adata = adata[:, common].copy()
adata.var["chromosome"]=("chr"+gp.loc[adata.var_names,"chrom"].astype(str)).values
adata.var["start"]=gp.loc[adata.var_names,"start"].astype(int).values
adata.var["end"]=gp.loc[adata.var_names,"end"].astype(int).values
# keep standard autosomes + X (exclude weird scaffolds)
ok = adata.var["chromosome"].isin(["chr"+str(i) for i in range(1,23)]+["chrX","chrY"])
adata = adata[:, ok].copy()
adata = adata[:, np.argsort(adata.var["chromosome"].values)].copy()  # helps ordering
print("genes with positions:", adata.n_vars)

lymph=["T","CTL","NK","B","ProB","Plasma"]
adata.obs["is_ref"]=(adata.obs["CellType"].isin(lymph) & (adata.obs["malignancy"]=="normal")).values

MALIG_MIN=100; REF_MIN=15
patients=[p for p in adata.obs["patient"].unique()
          if (adata.obs.patient.eq(p)&adata.obs.malignancy.eq("malignant")).sum()>=MALIG_MIN
          and (adata.obs.patient.eq(p)&adata.obs.is_ref).sum()>=REF_MIN]
print("analyzable patients:", len(patients), patients)

clone_rows=[]; summ_rows=[]
adata.obs["clone"]="NA"
for p in patients:
    sub = adata[(adata.obs.patient==p) & (adata.obs.malignancy.isin(["malignant","normal","unclear"]))].copy()
    sub.obs["cnv_ref"]=np.where(sub.obs["is_ref"],"ref","obs")
    cnv.tl.infercnv(sub, reference_key="cnv_ref", reference_cat=["ref"], window_size=100, step=10)
    cnv.tl.pca(sub, use_rep="cnv", key_added="cnv_pca", n_comps=min(30, sub.n_obs-1))
    cnv.pp.neighbors(sub, use_rep="cnv_pca")
    cnv.tl.leiden(sub, resolution=0.4)
    # per-cell CNV magnitude
    Xc = sub.obsm["X_cnv"]
    cell_score = np.asarray((Xc.multiply(Xc)).mean(axis=1)).ravel() if hasattr(Xc,"multiply") else (Xc**2).mean(axis=1)
    sub.obs["cnv_mag"]=cell_score
    ref_mag = sub.obs.loc[sub.obs.cnv_ref=="ref","cnv_mag"].mean()
    # per cluster
    cl = sub.obs.groupby("cnv_leiden")
    for g,idx in cl.groups.items():
        d=sub.obs.loc[idx]
        n=len(d); frac_malig=(d.malignancy=="malignant").mean()
        cmag=d.cnv_mag.mean(); ratio=cmag/ref_mag if ref_mag>0 else np.nan
        is_clone = (n>=30) and (frac_malig>=0.5) and (ratio>1.5)
        clone_rows.append(dict(patient=p,leiden=g,n=n,frac_malig=round(frac_malig,3),
                               cnv_ratio=round(ratio,3),is_malig_clone=is_clone))
        if is_clone:
            cells=d.index[d.malignancy=="malignant"]
            adata.obs.loc[cells,"clone"]=f"{p}_c{g}"
    # patient-level whole-malignant ratio
    malig_mag=sub.obs.loc[sub.obs.malignancy=="malignant","cnv_mag"].mean()
    nclones=sum(1 for r in clone_rows if r["patient"]==p and r["is_malig_clone"])
    summ_rows.append(dict(patient=p,n_malig=(sub.obs.malignancy=="malignant").sum(),
                          n_ref=(sub.obs.cnv_ref=="ref").sum(),
                          malig_cnv_ratio=round(malig_mag/ref_mag,3) if ref_mag>0 else np.nan,
                          n_clones=nclones))
    print(f"{p}: malig_ratio={malig_mag/ref_mag:.2f} clones={nclones}")

pd.DataFrame(clone_rows).to_csv("/tmp/aml/aml_clusters_all.csv",index=False)
pd.DataFrame(summ_rows).to_csv("/tmp/aml/aml_infercnv_summary.csv",index=False)
adata.write("/tmp/aml/aml_cnv.h5ad")
print("SUMMARY"); print(pd.DataFrame(summ_rows).to_string())
print("clones total:", (adata.obs.clone!="NA").sum(), "labels:", adata.obs.clone.nunique()-1)
