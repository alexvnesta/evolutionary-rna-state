import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial

a=ad.read_h5ad("/tmp/ddlps/ddlps_raw.h5ad")
gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"]); keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str); gp=gp.set_index("symbol")
for col in ["chromosome","start","end"]: a.var[col]=gp[col].reindex(a.var_names).values
a.layers["counts"]=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
a.layers["lognorm"]=a.X.copy()
genes=list(a.var_names)

markers={"immune_stromal_ref":["PTPRC","CD3D","LYZ","CD68","MS4A1","PECAM1","VWF"]}
pcounts=a.obs.patient.value_counts(); pats=[p for p in pcounts.index if pcounts[p]>=1500]

# accumulate pooled sums for malignant vs reference compartment
malig_sum=np.zeros(len(genes)); malig_n=0
ref_sum=np.zeros(len(genes)); ref_n=0
pb_rows=[]; pb_meta=[]
for pat in sorted(pats):
    ap=a[a.obs.patient==pat].copy(); ap.X=ap.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(ap,n_top_genes=2000)
    aph=ap[:,ap.var.highly_variable].copy(); sc.pp.scale(aph,max_value=10)
    sc.tl.pca(aph,n_comps=20); sc.pp.neighbors(aph,n_neighbors=15); sc.tl.leiden(aph,resolution=0.5,key_added="ct")
    ap.obs["ct"]=aph.obs["ct"].values
    refgenes=[g for g in markers["immune_stromal_ref"] if g in set(ap.var_names)]
    sc.tl.score_genes(ap,refgenes,score_name="ref_score")
    clust_ref=ap.obs.groupby("ct")["ref_score"].mean()
    ref_clusters=clust_ref[clust_ref>clust_ref.median()+0.1].index.tolist()
    if len(ref_clusters)<1: ref_clusters=[clust_ref.idxmax()]
    ap.obs["is_ref"]=ap.obs["ct"].isin(ref_clusters)
    n_ref=int(ap.obs["is_ref"].sum()); n_malig_cand=int((~ap.obs["is_ref"]).sum())
    if n_ref<50 or n_malig_cand<100: print(f"{pat}: skip",flush=True); continue
    ap.obs["ref_cat"]=np.where(ap.obs["is_ref"],"reference","malignant")
    cnv.tl.infercnv(ap,reference_key="ref_cat",reference_cat=["reference"],window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=20); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15); cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1)
    maligmask=(~ap.obs["is_ref"]).values; rc=ap.obs.loc[ap.obs.is_ref,"cnv_score"].mean()
    ct2=pd.crosstab(ap.obs["cnv_clone"],~ap.obs["is_ref"]); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct2.loc[c].get(True,0)>ct2.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    LN=ap.layers["lognorm"]
    # pooled compartment sums: malignant cells in called malignant clones; reference = is_ref cells
    malig_cells=ap.obs["cnv_clone"].isin(malig).values
    if malig_cells.sum()>0:
        malig_sum+=np.asarray(LN[malig_cells].sum(axis=0)).ravel(); malig_n+=int(malig_cells.sum())
    refmask=ap.obs["is_ref"].values
    if refmask.sum()>0:
        ref_sum+=np.asarray(LN[refmask].sum(axis=0)).ravel(); ref_n+=int(refmask.sum())
    for c in malig:
        m=(ap.obs.cnv_clone==c).values
        pb_rows.append(np.asarray(LN[m].mean(axis=0)).ravel())
        pb_meta.append({"patient":pat,"clone":f"{pat}_{c}","n_cells":int(m.sum())})
    print(f"{pat}: malig_clones={len(malig)} malig_cells={int(malig_cells.sum())} ref={int(refmask.sum())}",flush=True)

malig_mean=malig_sum/max(malig_n,1); ref_mean=ref_sum/max(ref_n,1)
comp=pd.DataFrame({"gene":genes,"malig_mean":malig_mean,"ref_mean":ref_mean,"spec":malig_mean-ref_mean})
comp.to_csv("/tmp/ddlps/ddlps_compartment_specificity.csv",index=False)
np.save("/tmp/ddlps/ddlps_clone_pb2.npy",np.vstack(pb_rows))
pd.DataFrame(pb_meta).to_csv("/tmp/ddlps/ddlps_clone_pb2_meta.csv",index=False)
print(f"malig_n={malig_n} ref_n={ref_n} specgenes>0.5={int((comp.spec>0.5).sum())} >1.0={int((comp.spec>1.0).sum())}",flush=True)
