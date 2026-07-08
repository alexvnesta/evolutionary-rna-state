import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial(fn,*it,**k): return [fn(*a) for a in zip(*it)]
_icnv.process_map=_serial

a=ad.read_h5ad("/tmp/luad/luad_raw.h5ad")
# gene positions
gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chr","start","end"]
gp=gp.dropna(subset=["symbol","chr","start"]); keep=[str(i) for i in range(1,23)]+["X","Y"]
gp=gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol").set_index("symbol")
gp["chromosome"]="chr"+gp["chr"].astype(str)
for c in ["chromosome","start","end"]: a.var[c]=gp[c].reindex(a.var_names).values
a.layers["counts"]=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); a.layers["lognorm"]=a.X.copy()
genes=list(a.var_names)

# LUAD hallmark programs
prog_defs={
 "epithelial_diff":["SFTPC","SFTPB","SFTPA1","NAPSA","SCGB1A1","MUC1","KRT8","KRT18"],
 "alveolar":["AGER","PDPN","SFTPC","HOPX"],
 "prolif":["MKI67","TOP2A","PCNA","CCNB1","CDK1","BIRC5"],
 "EMT_invasion":["VIM","FN1","SNAI2","ZEB1","CDH2","SPARC"],
 "interferon_MHC":["HLA-A","HLA-B","HLA-C","B2M","STAT1","IRF1","TAP1","GBP1","CXCL10"],
 "hypoxia_stress":["VEGFA","CA9","HIF1A","NDRG1","SLC2A1"],
}
ref_markers=["PTPRC","CD3D","CD3E","NKG7","LYZ","CD68","MS4A1","CD79A"]
pcounts=a.obs.patient.value_counts()
pats=[p for p in pcounts.index if pcounts[p]>=300]

summ=[]; val=[]; pb_rows=[]; pb_meta=[]; prog_rows=[]
for pat in sorted(pats):
    ap=a[a.obs.patient==pat].copy(); ap.X=ap.layers["lognorm"].copy()
    # reference = non-epithelial immune/stromal by author label
    ap.obs["is_ref"]=~ap.obs["cell_type"].isin(["Epithelial cells","Undetermined"])
    n_ref=int(ap.obs["is_ref"].sum()); n_epi=int((ap.obs["cell_type"]=="Epithelial cells").sum())
    if n_ref<50 or n_epi<100: print(f"{pat}: skip nref={n_ref} nepi={n_epi}",flush=True); continue
    ap.obs["ref_cat"]=np.where(ap.obs["is_ref"],"reference","malignant")
    cnv.tl.infercnv(ap,reference_key="ref_cat",reference_cat=["reference"],window_size=100,n_jobs=1)
    cnv.tl.pca(ap,n_comps=20); cnv.pp.neighbors(ap,use_rep="cnv_pca",n_neighbors=15); cnv.tl.leiden(ap,resolution=0.4,key_added="cnv_clone")
    Xc=ap.obsm["X_cnv"]; Xc=Xc.toarray() if hasattr(Xc,"toarray") else Xc
    ap.obs["cnv_score"]=np.abs(Xc).mean(1)
    rc=ap.obs.loc[ap.obs.is_ref,"cnv_score"].mean()
    epi_cnv=ap.obs.loc[ap.obs.cell_type=="Epithelial cells","cnv_score"].mean()
    ratio=epi_cnv/rc if rc>0 else np.nan
    ct2=pd.crosstab(ap.obs["cnv_clone"],ap.obs["cell_type"]=="Epithelial cells")
    es=ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    malig=[c for c in ap.obs["cnv_clone"].cat.categories
           if ct2.loc[c].get(True,0)>ct2.loc[c].sum()*0.5 and es[c]>1.5*rc and (ap.obs["cnv_clone"]==c).sum()>=30]
    summ.append({"patient":pat,"n_epi":n_epi,"n_ref":n_ref,"cnv_ratio":round(float(ratio),3),"n_malig_clones":len(malig)})
    print(f"{pat}: ratio={ratio:.2f} malig_clones={len(malig)} nepi={n_epi}",flush=True)
    LN=ap.layers["lognorm"]
    for c in malig:
        m=(ap.obs.cnv_clone==c).values
        pb_rows.append(np.asarray(LN[m].mean(axis=0)).ravel())
        pb_meta.append({"patient":pat,"clone":f"{pat}_{c}","n_cells":int(m.sum())})
    # validation: silhouette of clone partition in HVG-PCA space vs permutation null (multi-clone only)
    if len(malig)>=2:
        am=ap[ap.obs.cnv_clone.isin(malig)].copy(); am.X=am.layers["lognorm"].copy()
        sc.pp.highly_variable_genes(am,n_top_genes=2000); amh=am[:,am.var.highly_variable].copy()
        sc.pp.scale(amh,max_value=10); sc.tl.pca(amh,n_comps=20)
        from sklearn.metrics import silhouette_score
        lab=am.obs.cnv_clone.astype(str).values; emb=amh.obsm["X_pca"]
        sil=silhouette_score(emb,lab)
        null=[]; rng=np.random.default_rng(0)
        for _ in range(30): null.append(silhouette_score(emb,rng.permutation(lab)))
        null=np.array(null); z=(sil-null.mean())/(null.std()+1e-9)
        val.append({"patient":pat,"silhouette":round(float(sil),4),"null_mean":round(float(null.mean()),4),"z":round(float(z),2)})
        print(f"   val z={z:.1f} sil={sil:.3f}",flush=True)
    # programs per clone
    for c in malig:
        m=(ap.obs.cnv_clone==c).values
        row={"patient":pat,"clone":f"{pat}_{c}","n_cells":int(m.sum())}
        for pname,gl in prog_defs.items():
            gg=[g for g in gl if g in set(ap.var_names)]
            if gg: sc.tl.score_genes(ap,gg,score_name="_s"); row[pname]=float(ap.obs.loc[m,"_s"].mean())
            else: row[pname]=np.nan
        prog_rows.append(row)

pd.DataFrame(summ).to_csv("/tmp/luad/luad_infercnv_summary.csv",index=False)
pd.DataFrame(val).to_csv("/tmp/luad/luad_validation.csv",index=False)
np.save("/tmp/luad/luad_clone_pseudobulk.npy",np.vstack(pb_rows))
pd.DataFrame(pb_meta).to_csv("/tmp/luad/luad_clone_pseudobulk_meta.csv",index=False)
pd.Series(genes).to_csv("/tmp/luad/luad_pb_genes.csv",index=False)
pd.DataFrame(prog_rows).to_csv("/tmp/luad/luad_clone_programs.csv",index=False)
nmulti=sum(1 for s in summ if s["n_malig_clones"]>=2)
print(f"DONE: {len(summ)} patients, {nmulti} multi-clone, {len(pb_rows)} clones total",flush=True)
if val:
    zs=[v["z"] for v in val]; print(f"validation z: median={np.median(zs):.1f} range={min(zs):.1f}-{max(zs):.1f} above_null={sum(z>2 for z in zs)}/{len(zs)}",flush=True)
