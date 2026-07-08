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
for col in ["chromosome","start","end"]:
    a.var[col]=gp[col].reindex(a.var_names).values
a.X=a.layers["counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
lognorm=a.X.copy()
ref_types=["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]

# HVGs on malignant cells only, for a clean transcriptional-signal test
rows=[]
sil_real=[]; sil_null=[]
from sklearn.metrics import silhouette_score
for pat in sorted(a.obs["patient"].unique()):
    ap=a[a.obs["patient"]==pat].copy()
    refs=[t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    if not refs or (ap.obs["level1_celltype"]=="Epithelial").sum()<50: continue
    cnv.tl.infercnv(ap,reference_key="level1_celltype",reference_cat=refs,window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=min(20,ap.n_obs-1)); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15)
    cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1); epi=ap.obs["level1_celltype"]=="Epithelial"
    ct=pd.crosstab(ap.obs["cnv_clone"],epi); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    rc=ap.obs.loc[~epi,"cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct.loc[c].get(True,0)>ct.loc[c].sum()*0.5 and es[c]>1.5*rc]
    malig=[c for c in malig if (ap.obs["cnv_clone"]==c).sum()>=30]
    if len(malig)<2: 
        print(f"{pat}: <2 malig clones, skip validation"); continue
    # malignant cells, EXPRESSION space (HVG), NOT cnv space -> independent readout
    mcells=ap[ap.obs["cnv_clone"].isin(malig)].copy()
    sc.pp.highly_variable_genes(mcells,n_top_genes=2000)
    mm=mcells[:,mcells.var.highly_variable].copy()
    sc.pp.scale(mm,max_value=10); sc.tl.pca(mm,n_comps=20)
    lab=mm.obs["cnv_clone"].astype(str).values
    emb=mm.obsm["X_pca"]
    s_real=silhouette_score(emb,lab)
    # null: silhouette of RANDOM partition of same cluster sizes
    rngp=np.random.default_rng(42); nulls=[]
    for _ in range(50):
        perm=rngp.permutation(lab); nulls.append(silhouette_score(emb,perm))
    s_null=np.mean(nulls)
    rows.append({"patient":pat,"n_malig_clones":len(malig),"n_cells":mm.n_obs,
                 "silhouette_expr":s_real,"silhouette_null":s_null,"null_sd":np.std(nulls),
                 "z_over_null":(s_real-s_null)/(np.std(nulls)+1e-9)})
    print(f"{pat}: clones={len(malig)} cells={mm.n_obs} sil_expr={s_real:.3f} null={s_null:.3f}±{np.std(nulls):.3f} z={(s_real-s_null)/(np.std(nulls)+1e-9):.1f}")

pd.DataFrame(rows).to_csv("/tmp/cscc/phase0_validation.csv",index=False)
print("\nsaved phase0_validation.csv")
