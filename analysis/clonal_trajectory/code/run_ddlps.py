import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial
from sklearn.metrics import silhouette_score

a=ad.read_h5ad("/tmp/ddlps/ddlps_raw.h5ad")
# gene positions
gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"]); keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str); gp=gp.set_index("symbol")
for col in ["chromosome","start","end"]: a.var[col]=gp[col].reindex(a.var_names).values
a.layers["counts"]=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
a.layers["lognorm"]=a.X.copy()

# marker sets for compartment ID (immune/stromal = diploid reference; DDLPS malignant = adipocytic/mesenchymal)
markers={"Tcell":["CD3D","CD3E","CD2"],"Myeloid":["LYZ","CD68","AIF1","CD14"],
         "Bcell":["MS4A1","CD79A"],"Endo":["PECAM1","VWF","CDH5"],
         "immune_stromal_ref":["PTPRC","CD3D","LYZ","CD68","MS4A1","PECAM1","VWF"]}

val_rows=[]; summ_rows=[]
prog_defs={"adipocytic":["LPL","ADIPOQ","FABP4","PLIN1","CD36","PPARG"],
           "prolif":["MKI67","TOP2A","PCNA","CCNB1","CDK1"],
           "mesench_stem":["PDGFRA","PDGFRB","THY1","ENG","NT5E"],
           "MDM2_amp_region":["MDM2","CDK4","HMGA2","FRS2"],  # 12q amplicon hallmark of DDLPS
           "interferon_MHC":["HLA-A","HLA-B","HLA-C","B2M","STAT1","TAP1","CD74","HLA-DRA"],
           "EMT":["VIM","FN1","ZEB1","SNAI2","TWIST1"]}
prog_rows=[]

# choose patients with >=1500 cells
pcounts=a.obs.patient.value_counts()
pats=[p for p in pcounts.index if pcounts[p]>=1500]
print("patients analyzed:",len(pats))

for pat in sorted(pats):
    ap=a[a.obs.patient==pat].copy()
    ap.X=ap.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(ap,n_top_genes=2000)
    aph=ap[:,ap.var.highly_variable].copy(); sc.pp.scale(aph,max_value=10)
    sc.tl.pca(aph,n_comps=20); sc.pp.neighbors(aph,n_neighbors=15); sc.tl.leiden(aph,resolution=0.5,key_added="ct")
    ap.obs["ct"]=aph.obs["ct"].values
    # score immune/stromal reference markers per cluster
    refgenes=[g for g in markers["immune_stromal_ref"] if g in set(ap.var_names)]
    sc.tl.score_genes(ap,refgenes,score_name="ref_score")
    clust_ref=ap.obs.groupby("ct")["ref_score"].mean()
    # reference clusters = top ref_score (immune/stromal); malignant candidates = low ref_score
    ref_clusters=clust_ref[clust_ref>clust_ref.median()+0.1].index.tolist()
    if len(ref_clusters)<1:
        ref_clusters=[clust_ref.idxmax()]
    ap.obs["is_ref"]=ap.obs["ct"].isin(ref_clusters)
    n_ref=int(ap.obs["is_ref"].sum()); n_malig_cand=int((~ap.obs["is_ref"]).sum())
    if n_ref<50 or n_malig_cand<100:
        print(f"{pat}: skip (ref={n_ref}, malig_cand={n_malig_cand})"); continue
    ap.obs["ref_cat"]=np.where(ap.obs["is_ref"],"reference","malignant")
    cnv.tl.infercnv(ap,reference_key="ref_cat",reference_cat=["reference"],window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=20); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15); cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1)
    maligmask=~ap.obs["is_ref"]
    rc=ap.obs.loc[ap.obs.is_ref,"cnv_score"].mean()
    ct2=pd.crosstab(ap.obs["cnv_clone"],maligmask); es=ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories if ct2.loc[c].get(True,0)>ct2.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    ratio=ap.obs.loc[maligmask,"cnv_score"].mean()/rc
    summ_rows.append({"patient":pat,"n_cells":ap.n_obs,"n_ref":n_ref,"n_malig_clones":len(malig),
                      "n_malig_cells":int(ap.obs.cnv_clone.isin(malig).sum()),"cnv_ratio":ratio})
    # validation on multi-clone
    if len(malig)>=2:
        mc=ap[ap.obs.cnv_clone.isin(malig)].copy(); mc.X=mc.layers["lognorm"].copy()
        sc.pp.highly_variable_genes(mc,n_top_genes=2000); mm=mc[:,mc.var.highly_variable].copy()
        sc.pp.scale(mm,max_value=10); sc.tl.pca(mm,n_comps=20)
        lab=mm.obs.cnv_clone.astype(str).values; emb=mm.obsm["X_pca"]
        sr=silhouette_score(emb,lab); rng=np.random.default_rng(1)
        nulls=[silhouette_score(emb,rng.permutation(lab)) for _ in range(30)]
        val_rows.append({"patient":pat,"n_clones":len(malig),"sil":sr,"null":np.mean(nulls),"z":(sr-np.mean(nulls))/(np.std(nulls)+1e-9)})
        # programs per clone
        for prog,gl in prog_defs.items():
            gl2=[g for g in gl if g in set(ap.var_names)]
            if gl2: sc.tl.score_genes(ap,gl2,score_name=f"S_{prog}")
        for c in malig:
            m=(ap.obs.cnv_clone==c).values
            row={"patient":pat,"clone":f"{pat}_{c}","n_cells":int(m.sum())}
            for prog in prog_defs:
                if f"S_{prog}" in ap.obs: row[prog]=float(ap.obs.loc[m,f"S_{prog}"].mean())
            prog_rows.append(row)
    print(f"{pat}: cells={ap.n_obs} ref={n_ref} malig_clones={len(malig)} ratio={ratio:.2f}")

pd.DataFrame(summ_rows).to_csv("/tmp/ddlps/ddlps_infercnv_summary.csv",index=False)
pd.DataFrame(val_rows).to_csv("/tmp/ddlps/ddlps_validation.csv",index=False)
pd.DataFrame(prog_rows).to_csv("/tmp/ddlps/ddlps_clone_programs.csv",index=False)
print("\n=== DDLPS summary ===")
print(pd.DataFrame(summ_rows).round(2).to_string(index=False))
print("\nmulti-clone patients validated:",len(val_rows))
if val_rows: print(pd.DataFrame(val_rows).round(2).to_string(index=False))
