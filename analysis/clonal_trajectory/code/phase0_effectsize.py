import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial

a = ad.read_h5ad("/tmp/cscc/cscc_raw.h5ad")
gp = pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"])
keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str); gp=gp.set_index("symbol")
a.var["chromosome"]=gp["chromosome"].reindex(a.var_names).values
a.var["start"]=gp["start"].reindex(a.var_names).values
a.var["end"]=gp["end"].reindex(a.var_names).values

a.X=a.layers["counts"].copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
a.raw=a  # keep log-normalized expression for pseudobulk

ref_types=["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]

clone_pb=[]      # clone-level pseudobulk expression (log-norm mean over cells)
clone_meta=[]    # clone metadata
for pat in sorted(a.obs["patient"].unique()):
    ap=a[a.obs["patient"]==pat].copy()
    refs=[t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    n_epi=int((ap.obs["level1_celltype"]=="Epithelial").sum())
    if not refs or n_epi<50: 
        print(f"{pat}: skip"); continue
    cnv.tl.infercnv(ap,reference_key="level1_celltype",reference_cat=refs,window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=min(20,ap.n_obs-1))
    cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15)
    cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1)
    epi=ap.obs["level1_celltype"]=="Epithelial"
    ct=pd.crosstab(ap.obs["cnv_clone"],epi); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    rc=ap.obs.loc[~epi,"cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories
           if ct.loc[c].get(True,0)>ct.loc[c].sum()*0.5 and es[c]>1.5*rc]
    # pseudobulk each malignant clone (>=30 cells) from log-norm expression
    expr=ap.raw.X; expr=expr.toarray() if hasattr(expr,"toarray") else np.asarray(expr)
    for c in malig:
        mask=(ap.obs["cnv_clone"]==c).values
        if mask.sum()<30: continue
        clone_pb.append(expr[mask].mean(0))
        clone_meta.append({"patient":pat,"clone":f"{pat}_{c}","n_cells":int(mask.sum()),
                           "mean_cnv":float(ap.obs.loc[mask,"cnv_score"].mean())})
    print(f"{pat}: malig_clones_pb={sum(1 for m in clone_meta if m['patient']==pat)}")

PB=np.vstack(clone_pb)  # clones x genes (log-norm)
meta=pd.DataFrame(clone_meta)
genes=a.var_names.values
np.save("/tmp/cscc/clone_pseudobulk.npy",PB)
meta.to_csv("/tmp/cscc/clone_pseudobulk_meta.csv",index=False)
pd.Series(genes).to_csv("/tmp/cscc/pb_genes.csv",index=False)
print("\nclone pseudobulk matrix:",PB.shape,"| clones:",len(meta),"| patients:",meta.patient.nunique())
print(meta.groupby("patient").size().to_dict())
