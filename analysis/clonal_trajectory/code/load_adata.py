import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
import warnings; warnings.filterwarnings("ignore")

# metadata (per-cell, clean author annotations)
meta = pd.read_csv("GSE144236_patient_metadata_new.txt.gz", sep="\t", index_col=0)
print("meta:", meta.shape)
print("level1:", meta.level1_celltype.value_counts().to_dict())
print("patients:", meta.patient.value_counts().sort_index().to_dict())
print("tum.norm:", meta['tum.norm'].value_counts().to_dict())
print("level2 (top):", meta.level2_celltype.value_counts().head(12).to_dict())

# counts: rows=genes(+2 annotation rows), cols=cells. Read with pandas, drop annotation rows.
# First data rows are 'Patient' and 'Tissue...' — genes start after.
df = pd.read_csv("GSE144236_cSCC_counts.txt.gz", sep="\t", index_col=0)
print("\nraw counts frame:", df.shape, "| first idx:", list(df.index[:4]))
# drop the annotation rows if present
drop = [r for r in df.index if r in ("Patient","Tissue: 0=Normal, 1=Tumor")]
df = df.drop(index=drop)
print("after dropping annot rows:", df.shape)
# genes x cells -> transpose to cells x genes
X = df.T
# align to metadata cells
common = X.index.intersection(meta.index)
print("cells in counts:", X.shape[0], "| in meta:", meta.shape[0], "| common:", len(common))
X = X.loc[common]; m = meta.loc[common]
adata = ad.AnnData(X=X.values.astype(np.float32),
                   obs=m.copy(),
                   var=pd.DataFrame(index=X.columns))
adata.layers["counts"] = adata.X.copy()
adata.write("/tmp/cscc/cscc_raw.h5ad")
print("\nsaved AnnData:", adata.shape)
print("obs cols:", list(adata.obs.columns))
