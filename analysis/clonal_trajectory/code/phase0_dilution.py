import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial

a = ad.read_h5ad("/tmp/cscc/cscc_raw.h5ad")
gp = pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"]); keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str); gp=gp.set_index("symbol")
for col in ["chromosome","start","end"]: a.var[col]=gp[col].reindex(a.var_names).values
a.X=a.layers["counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
lognorm=a.X.copy()
a.X=a.layers["counts"].copy()  # keep raw counts for realistic remixing
ref_types=["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]

# get malignant clone labels (rerun infercnv, store labels per cell)
labels = pd.Series(index=a.obs_names, dtype=object)
for pat in sorted(a.obs["patient"].unique()):
    ap=a[a.obs["patient"]==pat].copy()
    ap.X=lognorm[(a.obs["patient"]==pat).values]  # infercnv on lognorm
    refs=[t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    if not refs or (ap.obs["level1_celltype"]=="Epithelial").sum()<50: continue
    cnv.tl.infercnv(ap,reference_key="level1_celltype",reference_cat=refs,window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=min(20,ap.n_obs-1)); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15)
    cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1); epi=ap.obs["level1_celltype"]=="Epithelial"
    ct=pd.crosstab(ap.obs["cnv_clone"],epi); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean(); rc=ap.obs.loc[~epi,"cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct.loc[c].get(True,0)>ct.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    for c in malig:
        cells=ap.obs_names[ap.obs["cnv_clone"]==c]
        labels[cells]=f"{pat}_{c}"

a.obs["clone"]=labels
raw=a.layers["counts"]

# DILUTION TEST: for each patient, build "bulk" = clone malignant cells + patient's own non-malignant cells
# at target purity p. Measure whether clone identity survives dilution (silhouette in expr space).
from sklearn.metrics import silhouette_score
rng=np.random.default_rng(0)
def make_bulk(cell_idx, n_reps=40, n_cells=200):
    # sample n_cells with replacement, sum raw counts, CPM-log normalize -> one pseudobulk replicate
    reps=[]
    for _ in range(n_reps):
        s=rng.choice(cell_idx,size=n_cells,replace=True)
        v=np.asarray(raw[s].sum(0)).ravel()
        v=np.log1p(v/ v.sum()*1e4)
        reps.append(v)
    return np.array(reps)

purities=[1.0,0.6,0.3,0.1]
rows=[]
for pat in sorted(a.obs["patient"].dropna().unique()):
    pmask=(a.obs["patient"]==pat).values
    clones=[c for c in a.obs.loc[pmask,"clone"].dropna().unique()]
    if len(clones)<2: continue
    nonmalig=np.where(pmask & a.obs["clone"].isna().values & ~a.obs["level1_celltype"].isin(["Multiplet"]).values)[0]
    if len(nonmalig)<50: 
        # fall back to any non-clone cell
        nonmalig=np.where(pmask & a.obs["clone"].isna().values)[0]
    for pur in purities:
        X=[]; y=[]
        for c in clones:
            cidx=np.where(pmask & (a.obs["clone"]==c).values)[0]
            if len(cidx)<30: continue
            for _ in range(30):
                nm=int(200*(1-pur)); nc=200-nm
                s_c=rng.choice(cidx,nc,replace=True)
                s_n=rng.choice(nonmalig,max(nm,0),replace=True) if nm>0 else np.array([],dtype=int)
                s=np.concatenate([s_c,s_n]).astype(int)
                v=np.asarray(raw[s].sum(0)).ravel(); v=np.log1p(v/v.sum()*1e4)
                X.append(v); y.append(c)
        X=np.array(X); y=np.array(y)
        # HVG-ish: top variance genes
        vv=X.var(0); topg=np.argsort(vv)[-2000:]
        Xs=X[:,topg]; Xs=(Xs-Xs.mean(0))/(Xs.std(0)+1e-9)
        from sklearn.decomposition import PCA
        emb=PCA(n_components=min(20,Xs.shape[0]-1)).fit_transform(Xs)
        sil=silhouette_score(emb,y)
        rows.append({"patient":pat,"purity":pur,"n_clones":len(set(y)),"silhouette":sil})
    print(pat, "done")

df=pd.DataFrame(rows); df.to_csv("/tmp/cscc/phase0_dilution.csv",index=False)
print("\n=== Class separability vs tumor purity ===")
print(df.groupby("purity")["silhouette"].agg(["mean","std","count"]).round(3))
