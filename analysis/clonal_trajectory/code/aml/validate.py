import os, warnings, numpy as np, pandas as pd, scanpy as sc
np.random.seed(0); warnings.filterwarnings("ignore")
from sklearn.metrics import silhouette_score

adata=sc.read_h5ad("/tmp/aml/aml_cnv.h5ad")  # has clone labels; X is raw counts
summ=pd.read_csv("/tmp/aml/aml_infercnv_summary.csv")
multiclone=summ[summ.n_clones>=2].patient.tolist()
print("multiclone patients:",multiclone)

# ---- normalize expression once (log1p CPM) for HVG-space validation & pseudobulk ----
adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
adata.raw=adata

N_PERM=50
val_rows=[]
for p in multiclone:
    mal=adata[(adata.obs.patient==p)&(adata.obs.clone!="NA")].copy()
    labels=mal.obs["clone"].astype(str).values
    if len(set(labels))<2: continue
    # independent expression space: HVG PCA (NOT cnv features)
    sc.pp.highly_variable_genes(mal,n_top_genes=2000)
    m2=mal[:,mal.var.highly_variable].copy()
    sc.pp.scale(m2,max_value=10)
    sc.tl.pca(m2,n_comps=min(30,m2.n_obs-1))
    X=m2.obsm["X_pca"]
    obs_sil=silhouette_score(X,labels)
    # within-patient permutation null (shuffle clone labels)
    rng=np.random.default_rng(0); null=[]
    for _ in range(N_PERM):
        null.append(silhouette_score(X,rng.permutation(labels)))
    null=np.array(null)
    z=(obs_sil-null.mean())/null.std() if null.std()>0 else np.nan
    val_rows.append(dict(patient=p,n_clones=len(set(labels)),n_cells=mal.n_obs,
                         silhouette=round(obs_sil,4),null_mean=round(null.mean(),4),
                         null_sd=round(null.std(),4),z=round(z,3)))
    print(f"{p}: sil={obs_sil:.3f} null={null.mean():.3f}+-{null.std():.3f} z={z:.2f}")

vdf=pd.DataFrame(val_rows); vdf.to_csv("/tmp/aml/aml_validation.csv",index=False)
print(vdf.to_string())
print("median z:",round(vdf.z.median(),3),"n z>2:",int((vdf.z>2).sum()),"/",len(vdf))
