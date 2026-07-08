import os, warnings, numpy as np, pandas as pd, scanpy as sc
np.random.seed(0); warnings.filterwarnings("ignore")
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import pearsonr

adata=sc.read_h5ad("/tmp/aml/aml_cnv.h5ad")
adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)

clones=sorted([c for c in adata.obs.clone.unique() if c!="NA"])
print("n clones:",len(clones))

# ---- pseudobulk each clone (mean log-norm) ----
genes=adata.var_names.to_numpy()
pb=np.zeros((len(clones),adata.n_vars),dtype=np.float32)
meta=[]
for i,c in enumerate(clones):
    m=adata.obs.clone==c
    pb[i]=np.asarray(adata[m].X.mean(0)).ravel()
    meta.append(dict(clone=c,patient=c.split("_c")[0],n=int(m.sum())))
meta=pd.DataFrame(meta)

# ---- effect size: within vs between patient separation in top-variance pseudobulk space ----
v=pb.var(0); top=np.argsort(v)[-2000:]
Z=(pb[:,top]-pb[:,top].mean(0))/(pb[:,top].std(0)+1e-9)
from scipy.spatial.distance import pdist, squareform
D=squareform(pdist(Z))  # euclidean per-gene => divide by sqrt(ngenes) for per-gene SD
D_pg=D/np.sqrt(Z.shape[1])
pat=meta.patient.values
wi=[]; bt=[]
for i in range(len(clones)):
    for j in range(i+1,len(clones)):
        (wi if pat[i]==pat[j] else bt).append(D_pg[i,j])
print(f"within-patient subclone sep: {np.mean(wi):.3f} SD/gene (n={len(wi)})")
print(f"between-patient sep: {np.mean(bt):.3f} SD/gene (n={len(bt)})")

# ---- program scores per clone ----
progs={
 "interferon_MHC":["HLA-A","HLA-B","HLA-C","B2M","STAT1","TAP1","CD74","HLA-DRA"],
 "HSC_progenitor":["CD34","KIT","MEIS1","HOXA9"],
 "myeloid_diff":["CD14","LYZ","MPO","ITGAM","CEBPE"],
 "cell_cycle":["MKI67","TOP2A","PCNA"],
}
gidx={g:i for i,g in enumerate(genes)}
def score(genelist):
    ii=[gidx[g] for g in genelist if g in gidx]
    return pb[:,ii].mean(1)
prog_mat=pd.DataFrame({k:score(v) for k,v in progs.items()})
for k,v in progs.items():
    miss=[g for g in v if g not in gidx]
    if miss: print("missing from",k,":",miss)

# ---- recurrent classes: Ward on programs, k=2 ----
S=(prog_mat-prog_mat.mean())/prog_mat.std()
L=linkage(S.values,method="ward")
cls=fcluster(L,2,criterion="maxclust")
meta["class_k2"]=cls
res=pd.concat([meta.reset_index(drop=True),prog_mat.reset_index(drop=True)],axis=1)
res.to_csv("/tmp/aml/aml_clone_programs_classified.csv",index=False)

# ---- immune coupling: correlate each program with IFN/MHC across clones ----
ifn=prog_mat["interferon_MHC"].values
corr_rows=[]
for k in ["HSC_progenitor","myeloid_diff","cell_cycle"]:
    r,pv=pearsonr(prog_mat[k].values,ifn)
    corr_rows.append(dict(program=k,r_vs_IFN_MHC=round(r,3),p=round(pv,4)))
    print(f"{k} vs IFN/MHC: r={r:.3f} p={pv:.4f}")
# class-level IFN/MHC difference
for cl in [1,2]:
    print(f"class {cl}: n={int((cls==cl).sum())} IFN/MHC mean={ifn[cls==cl].mean():.3f}")

# save pseudobulk artifacts
np.save("/tmp/aml/aml_clone_pseudobulk.npy",pb)
pd.Series(genes).to_csv("/tmp/aml/aml_pb_genes.csv",index=False,header=["gene"])
meta_out=meta[["clone","patient","class_k2"]].copy()
meta_out.to_csv("/tmp/aml/aml_clone_pseudobulk_meta.csv",index=False)
pd.DataFrame(corr_rows).to_csv("/tmp/aml/aml_immune_coupling.csv",index=False)

# effect size table
pd.DataFrame([dict(within_patient_SD=round(np.mean(wi),4),
                   between_patient_SD=round(np.mean(bt),4),
                   n_within=len(wi),n_between=len(bt))]).to_csv("/tmp/aml/aml_effectsize.csv",index=False)
print("class counts:",dict(pd.Series(cls).value_counts()))
print("saved all program artifacts")
