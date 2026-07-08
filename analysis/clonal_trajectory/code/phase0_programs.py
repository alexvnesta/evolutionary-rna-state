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
for col in ["chromosome","start","end"]: a.var[col]=gp[col].reindex(a.var_names).values
a.X=a.layers["counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
ref_types=["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]

# hallmark-ish marker sets relevant to cSCC clonal biology + immune phenotype
progs = {
 "epithelial_diff":["KRT1","KRT10","SBSN","CALML5","KRTDAP","FLG","LOR","CDKN2A"],
 "basal_stem":["KRT5","KRT14","COL17A1","TP63","ITGA6","ITGB1","KRT15"],
 "EMT_invasion":["VIM","FN1","MMP1","MMP10","LAMC2","ITGA5","SNAI2","TGFBI"],
 "cell_cycle":["MKI67","TOP2A","PCNA","CCNB1","CDK1","CENPF","UBE2C"],
 "interferon_MHC":["HLA-A","HLA-B","HLA-C","B2M","STAT1","IRF1","TAP1","HLA-DRA","CD74"],
 "hypoxia_stress":["VEGFA","HILPDA","NDRG1","BNIP3","LDHA","SLC2A1"],
}
allmk=sorted({g for v in progs.values() for g in v if g in set(a.var_names)})

rows=[]
for pat in sorted(a.obs["patient"].unique()):
    ap=a[a.obs["patient"]==pat].copy()
    refs=[t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    if not refs or (ap.obs["level1_celltype"]=="Epithelial").sum()<50: continue
    cnv.tl.infercnv(ap,reference_key="level1_celltype",reference_cat=refs,window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=min(20,ap.n_obs-1)); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15)
    cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1); epi=ap.obs["level1_celltype"]=="Epithelial"
    ct=pd.crosstab(ap.obs["cnv_clone"],epi); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean(); rc=ap.obs.loc[~epi,"cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct.loc[c].get(True,0)>ct.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    # score programs per clone
    for prog,gl in progs.items():
        gl2=[g for g in gl if g in set(ap.var_names)]
        sc.tl.score_genes(ap,gl2,score_name=f"S_{prog}")
    for c in malig:
        m=(ap.obs["cnv_clone"]==c).values
        row={"patient":pat,"clone":f"{pat}_{c}","n_cells":int(m.sum()),"cnv":float(ap.obs.loc[m,"cnv_score"].mean())}
        for prog in progs: row[prog]=float(ap.obs.loc[m,f"S_{prog}"].mean())
        rows.append(row)
    print(f"{pat}: scored {len(malig)} clones")

df=pd.DataFrame(rows)
df.to_csv("/tmp/cscc/clone_programs.csv",index=False)
print("\nsaved clone_programs.csv:",df.shape)
print(df.round(3).to_string(index=False))
