import os, glob
import numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
import anndata as ad

RAW="/tmp/aml/raw"
np.random.seed(0)

# map sample -> anno file (GSM differs between dem/anno)
anno_map={}
for f in glob.glob(os.path.join(RAW,"*.anno.txt.gz")):
    s = os.path.basename(f).split("_",1)[1].replace(".anno.txt.gz","")
    anno_map[s]=f

dem_files = sorted(glob.glob(os.path.join(RAW,"*-D0.dem.txt.gz")))
print("n D0 dem files:", len(dem_files))
adatas=[]
for dem in dem_files:
    rest = os.path.basename(dem).split("_",1)[1]
    sample = rest.replace(".dem.txt.gz","")
    patient = sample.split("-")[0]
    df = pd.read_csv(dem, sep="\t", index_col=0)
    X = sparse.csr_matrix(df.T.values.astype(np.float32))
    cells = df.columns.to_numpy(); genes = df.index.to_numpy()
    anno = pd.read_csv(anno_map[sample], sep="\t", index_col=0).reindex(cells)
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    a.obs["patient"]=patient; a.obs["sample"]=sample
    a.obs["CellType"]=anno["CellType"].astype(str).values
    a.obs["malignancy"]=anno["PredictionRefined"].astype(str).values
    a.obs["CyclingScore"]=pd.to_numeric(anno["CyclingScore"],errors="coerce").values
    a.obs["MutTranscripts"]=anno["MutTranscripts"].astype(str).values
    a.obs["WtTranscripts"]=anno["WtTranscripts"].astype(str).values
    adatas.append(a)
    print(f"{sample}: {a.n_obs} x {a.n_vars}")

adata = ad.concat(adatas, join="outer", index_unique=None, fill_value=0)
adata.obs_names_make_unique()
print("concat:", adata.shape)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print("after QC:", adata.shape)
print(adata.obs["malignancy"].value_counts().to_dict())
print("patients:", adata.obs["patient"].nunique())
adata.write("/tmp/aml/aml_raw.h5ad")
print("saved")
