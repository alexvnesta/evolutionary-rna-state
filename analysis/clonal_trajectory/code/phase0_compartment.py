import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

a=ad.read_h5ad("/tmp/cscc/cscc_raw.h5ad")
gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"]); keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str); gp=gp.set_index("symbol")
for col in ["chromosome","start","end"]: a.var[col]=gp[col].reindex(a.var_names).values
a.X=a.layers["counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); a.layers["lognorm"]=a.X.copy()
ref_types=["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]

# get clone labels + define malignant-compartment-specific genes (high in Epithelial, low in immune/stromal)
labels=pd.Series(index=a.obs_names,dtype=object)
epi_mask_all=(a.obs["level1_celltype"]=="Epithelial").values
nonepi_mask_all=a.obs["level1_celltype"].isin(ref_types).values
# compartment specificity score per gene: mean in epithelial minus mean in reference
Xln=a.layers["lognorm"]
epi_mean=np.asarray(Xln[epi_mask_all].mean(0)).ravel()
ref_mean=np.asarray(Xln[nonepi_mask_all].mean(0)).ravel()
malig_specificity=epi_mean-ref_mean  # positive = malignant-compartment-enriched
a.var["malig_spec"]=malig_specificity

for pat in sorted(a.obs["patient"].unique()):
    ap=a[a.obs["patient"]==pat].copy(); ap.X=ap.layers["lognorm"].copy()
    refs=[t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    if not refs or (ap.obs["level1_celltype"]=="Epithelial").sum()<50: continue
    cnv.tl.infercnv(ap,reference_key="level1_celltype",reference_cat=refs,window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=min(20,ap.n_obs-1)); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15); cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1); epi=ap.obs["level1_celltype"]=="Epithelial"
    ct=pd.crosstab(ap.obs["cnv_clone"],epi); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean(); rc=ap.obs.loc[~epi,"cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct.loc[c].get(True,0)>ct.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    for c in malig:
        labels[ap.obs_names[ap.obs["cnv_clone"]==c]]=f"{pat}_{c}"
a.obs["clone"]=labels

# dilution test: ALL top-variance genes vs MALIGNANT-SPECIFIC genes
raw=a.layers["counts"]; rng=np.random.default_rng(0)
purities=[1.0,0.6,0.3,0.1]
# malignant-specific gene set: top 2000 by specificity (and expressed)
malig_genes=np.where((a.var["malig_spec"]>0.5))[0]
print(f"malignant-compartment-specific genes (spec>0.5): {len(malig_genes)}")

rows=[]
for pat in sorted(a.obs["patient"].dropna().unique()):
    pmask=(a.obs["patient"]==pat).values
    clones=[c for c in a.obs.loc[pmask,"clone"].dropna().unique()]
    if len(clones)<2: continue
    nonmalig=np.where(pmask & a.obs["clone"].isna().values)[0]
    for pur in purities:
        X=[]; y=[]
        for c in clones:
            cidx=np.where(pmask & (a.obs["clone"]==c).values)[0]
            if len(cidx)<30: continue
            for _ in range(30):
                nm=int(200*(1-pur)); nc=200-nm
                s=np.concatenate([rng.choice(cidx,nc,replace=True),
                                  rng.choice(nonmalig,max(nm,0),replace=True) if nm>0 else np.array([],dtype=int)]).astype(int)
                v=np.asarray(raw[s].sum(0)).ravel(); v=np.log1p(v/v.sum()*1e4)
                X.append(v); y.append(c)
        X=np.array(X); y=np.array(y)
        def sil_on(genes):
            Xs=X[:,genes]; vv=Xs.var(0); tp=np.argsort(vv)[-2000:]
            Xs=Xs[:,tp]; Xs=(Xs-Xs.mean(0))/(Xs.std(0)+1e-9)
            emb=PCA(n_components=min(20,Xs.shape[0]-1)).fit_transform(Xs)
            return silhouette_score(emb,y)
        all_genes=np.arange(X.shape[1])
        rows.append({"patient":pat,"purity":pur,"sil_allgenes":sil_on(all_genes),"sil_maligspec":sil_on(malig_genes)})
    print(pat,"done")
df=pd.DataFrame(rows); df.to_csv("/tmp/cscc/phase0_compartment.csv",index=False)
print("\n=== dilution: all-genes vs malignant-specific-genes ===")
print(df.groupby("purity")[["sil_allgenes","sil_maligspec"]].mean().round(3))
