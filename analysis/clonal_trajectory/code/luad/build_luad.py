import gzip, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse
import warnings, time; warnings.filterwarnings("ignore"); t0=time.time()

anno=pd.read_csv("/tmp/luad/cell_annotation.txt.gz",sep="\t").set_index("Index")
tumor_origins=["tLung","mBrain","mLN","PE","tL/B"]
epi=anno[(anno["Cell_type"]=="Epithelial cells")&(anno["Sample_Origin"].isin(tumor_origins))]
good=epi.groupby("Sample").size(); good=good[good>=200].index.tolist()
keep_set=set(anno[anno["Sample"].isin(good)].index)
print(f"samples={len(good)} target_cells={len(keep_set)}",flush=True)

# stream: first line = cell barcodes; each subsequent line = gene + counts
f=gzip.open("/tmp/luad/raw_UMI.txt.gz","rt")
header=f.readline().rstrip("\n").split("\t")
cells=header[1:]
keep_idx=np.array([i for i,c in enumerate(cells) if c in keep_set])
kept_cells=[cells[i] for i in keep_idx]
print(f"kept cols={len(keep_idx)}",flush=True)

genes=[]; data=[]; rows=[]; cols=[]; gi=0
for line in f:
    parts=line.rstrip("\n").split("\t")
    genes.append(parts[0])
    vals=np.array(parts[1:],dtype=np.float32)[keep_idx]
    nz=np.nonzero(vals)[0]
    if len(nz):
        data.append(vals[nz]); cols.append(nz); rows.append(np.full(len(nz),gi))
    gi+=1
    if gi%5000==0: print(f"  {gi} genes, {time.time()-t0:.0f}s",flush=True)
f.close()
data=np.concatenate(data); rows=np.concatenate(rows); cols=np.concatenate(cols)
M=sparse.csr_matrix((data,(rows,cols)),shape=(len(genes),len(keep_idx)))  # genes x cells
print(f"sparse {M.shape} nnz={M.nnz} in {time.time()-t0:.0f}s",flush=True)
a=ad.AnnData(X=M.T.tocsr())  # cells x genes
a.obs_names=kept_cells; a.var_names=genes
meta=anno.loc[a.obs_names]
for c,src in [("patient","Sample"),("origin","Sample_Origin"),("cell_type","Cell_type"),("subtype","Cell_subtype")]:
    a.obs[c]=meta[src].values
a.write("/tmp/luad/luad_raw.h5ad")
print(f"saved {a.shape} patients={a.obs.patient.nunique()} in {time.time()-t0:.0f}s",flush=True)
print(a.obs.cell_type.value_counts().to_string(),flush=True)
