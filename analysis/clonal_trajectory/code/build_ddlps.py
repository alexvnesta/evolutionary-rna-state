import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix
import warnings; warnings.filterwarnings("ignore")

print("reading mtx...")
M = mmread("/tmp/ddlps/GSE221493_sc_LPS_atlas_raw_counts_matrix.mtx.gz").tocsr()  # genes x cells
print("mtx:", M.shape)
barc = pd.read_csv("/tmp/ddlps/GSE221493_sc_LPS_atlas_raw_counts_barcodes.tsv.gz",header=None)[0].values
feat = pd.read_csv("/tmp/ddlps/GSE221493_sc_LPS_atlas_raw_counts_features.tsv.gz",header=None,sep="\t")[0].values
print("barcodes:",len(barc),"features:",len(feat))

# genes x cells -> cells x genes
X = csr_matrix(M.T)
patient = pd.Series(barc).str.extract(r'-(\d+)$')[0].values
adata = ad.AnnData(X=X, obs=pd.DataFrame({"barcode":barc,"patient":["P"+p for p in patient]},index=barc),
                   var=pd.DataFrame(index=feat))
adata.var_names_make_unique()
adata.layers["counts"]=adata.X.copy()
# basic QC
sc.pp.filter_cells(adata,min_genes=200); sc.pp.filter_genes(adata,min_cells=10)
adata.obs["n_counts"]=np.asarray(adata.layers["counts"].sum(1)).ravel() if False else np.asarray(adata.X.sum(1)).ravel()
print("after QC:",adata.shape)
print("patients:",pd.Series(adata.obs.patient).value_counts().sort_index().to_dict())
adata.write("/tmp/ddlps/ddlps_raw.h5ad")
print("saved ddlps_raw.h5ad")
