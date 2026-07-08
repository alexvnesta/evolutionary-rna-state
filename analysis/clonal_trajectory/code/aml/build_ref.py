import os, glob, numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
import anndata as ad
np.random.seed(0)
RAW="/tmp/aml/raw"
anno_map={os.path.basename(f).split("_",1)[1].replace(".anno.txt.gz",""):f
          for f in glob.glob(os.path.join(RAW,"*.anno.txt.gz"))}
# BM samples -> all cells are diploid reference
bm_dems=sorted(glob.glob(os.path.join(RAW,"*_BM*.dem.txt.gz")))
print("BM dem files:", [os.path.basename(x) for x in bm_dems])
refs=[]
for dem in bm_dems:
    s=os.path.basename(dem).split("_",1)[1].replace(".dem.txt.gz","")
    df=pd.read_csv(dem,sep="\t",index_col=0)
    X=sparse.csr_matrix(df.T.values.astype(np.float32))
    a=ad.AnnData(X=X,obs=pd.DataFrame(index=df.columns.to_numpy()),
                 var=pd.DataFrame(index=df.index.to_numpy()))
    a.obs["refsample"]=s
    refs.append(a); print(f"{s}: {a.n_obs} cells")
bm=ad.concat(refs,join="outer",fill_value=0); bm.obs_names_make_unique()
sc.pp.filter_cells(bm,min_genes=200)
print("BM pooled ref:", bm.shape)
bm.write("/tmp/aml/bm_ref.h5ad")
