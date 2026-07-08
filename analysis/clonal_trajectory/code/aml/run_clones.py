import os, warnings, numpy as np, pandas as pd, scanpy as sc
np.random.seed(0); warnings.filterwarnings("ignore")
import anndata as ad, infercnvpy as cnv
import infercnvpy.tl._infercnv as _icnv
def _serial_map(fn,*its,**kw): return [fn(*a) for a in zip(*its)]
_icnv.process_map=_serial_map

adata=sc.read_h5ad("/tmp/aml/aml_raw.h5ad")
bm=sc.read_h5ad("/tmp/aml/bm_ref.h5ad")
bm=bm[bm.obs.refsample.isin(["BM1","BM2","BM3","BM4"])].copy()
common=adata.var_names.intersection(bm.var_names)
adata=adata[:,common].copy(); bm=bm[:,common].copy()

gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chrom","start","end"]
gp=gp.drop_duplicates("symbol").set_index("symbol")
def prep(A):
    c=A.var_names.intersection(gp.index); A=A[:,c].copy()
    A.var["chromosome"]=("chr"+gp.loc[A.var_names,"chrom"].astype(str)).values
    A.var["start"]=gp.loc[A.var_names,"start"].astype(int).values
    A.var["end"]=gp.loc[A.var_names,"end"].astype(int).values
    ok=A.var["chromosome"].isin(["chr"+str(i) for i in range(1,23)]+["chrX"])
    return A[:,ok].copy()

lymph=["T","CTL","NK","B","ProB","Plasma"]
lymph_ref=adata[(adata.obs.CellType.isin(lymph))&(adata.obs.malignancy=="normal")].copy()

MALIG_MIN=100; CLONE_MIN=30
patients=[p for p in adata.obs.patient.unique()
          if (adata.obs.patient.eq(p)&adata.obs.malignancy.eq("malignant")).sum()>=MALIG_MIN]

clone_rows=[]; summ_rows=[]
adata.obs["clone"]="NA"
store={}   # patient -> (cnv_pca of malignant, obs_names, expr_hvg_pca)
for p in patients:
    mal=adata[(adata.obs.patient==p)&(adata.obs.malignancy=="malignant")].copy()
    mal.obs["role"]="obs"
    r1=bm.copy(); r1.obs["role"]="ref"
    r2=lymph_ref[lymph_ref.obs.patient!=p].copy(); r2.obs["role"]="ref"
    sub=ad.concat([mal,r1,r2],join="inner"); sub=prep(sub)
    sc.pp.normalize_total(sub,target_sum=1e4); sc.pp.log1p(sub)
    cnv.tl.infercnv(sub,reference_key="role",reference_cat=["ref"],window_size=100,step=10)
    Xc=sub.obsm["X_cnv"]; mag=np.asarray((Xc.multiply(Xc)).mean(1)).ravel()
    sub.obs["cnv_mag"]=mag
    ref_mag=sub.obs.loc[sub.obs.role=="ref","cnv_mag"].mean()
    obs=sub[sub.obs.role=="obs"].copy()
    cnv.tl.pca(obs,use_rep="cnv",key_added="cnv_pca",n_comps=min(30,obs.n_obs-1))
    cnv.pp.neighbors(obs,use_rep="cnv_pca")
    cnv.tl.leiden(obs,resolution=0.4)
    obs.obs["cnv_mag"]=sub.obs.loc[obs.obs_names,"cnv_mag"].values
    malig_ratio=obs.obs["cnv_mag"].mean()/ref_mag
    # subclones = leiden clusters >=30 cells
    vc=obs.obs["cnv_leiden"].value_counts()
    keep=vc[vc>=CLONE_MIN].index.tolist()
    nclones=len(keep)
    for g in vc.index:
        n=int(vc[g]); ratio=obs.obs.loc[obs.obs.cnv_leiden==g,"cnv_mag"].mean()/ref_mag
        iscl=(n>=CLONE_MIN)
        clone_rows.append(dict(patient=p,leiden=str(g),n=n,cnv_ratio=round(ratio,3),is_clone=iscl))
        if iscl:
            cells=obs.obs_names[obs.obs.cnv_leiden==g]
            adata.obs.loc[cells,"clone"]=f"{p}_c{g}"
    summ_rows.append(dict(patient=p,n_malig=obs.n_obs,n_ref=int((sub.obs.role=="ref").sum()),
                          malig_cnv_ratio=round(malig_ratio,3),n_clones=nclones,
                          aneuploid=bool(malig_ratio>1.5)))
    # save cnv_pca for the multiclone validation
    store[p]=dict(names=obs.obs_names.to_numpy(),
                  cnv_pca=obs.obsm["X_cnv_pca"],
                  cnv_leiden=obs.obs["cnv_leiden"].to_numpy())
    print(f"{p}: ratio={malig_ratio:.2f} nclones={nclones} sizes={vc.to_dict()}")

pd.DataFrame(clone_rows).to_csv("/tmp/aml/aml_clusters_all.csv",index=False)
sdf=pd.DataFrame(summ_rows); sdf.to_csv("/tmp/aml/aml_infercnv_summary.csv",index=False)
np.save("/tmp/aml/cnv_store.npy",store,allow_pickle=True)
adata.write("/tmp/aml/aml_cnv.h5ad")
print(sdf.to_string())
print("multiclone patients:",int((sdf.n_clones>=2).sum()),
      "aneuploid:",int(sdf.aneuploid.sum()),
      "total clones:",adata.obs.clone.nunique()-1)
