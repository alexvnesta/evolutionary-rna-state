import os, warnings, numpy as np, pandas as pd, scanpy as sc
np.random.seed(0); warnings.filterwarnings("ignore")
import anndata as ad, infercnvpy as cnv
import infercnvpy.tl._infercnv as _icnv
def _serial_map(fn,*its,**kw): return [fn(*a) for a in zip(*its)]
_icnv.process_map=_serial_map

adata=sc.read_h5ad("/tmp/aml/aml_raw.h5ad")
bm=sc.read_h5ad("/tmp/aml/bm_ref.h5ad")
# unsorted normal marrow only (BM1-4), exclude CD34-sorted BM5
bm=bm[bm.obs.refsample.isin(["BM1","BM2","BM3","BM4"])].copy()

# gene positions -> restrict + annotate
gp=pd.read_csv("/tmp/cscc/gene_pos.csv"); gp.columns=["symbol","chrom","start","end"]
gp=gp.drop_duplicates("symbol").set_index("symbol")
def prep(A):
    c=A.var_names.intersection(gp.index); A=A[:,c].copy()
    A.var["chromosome"]=("chr"+gp.loc[A.var_names,"chrom"].astype(str)).values
    A.var["start"]=gp.loc[A.var_names,"start"].astype(int).values
    A.var["end"]=gp.loc[A.var_names,"end"].astype(int).values
    ok=A.var["chromosome"].isin(["chr"+str(i) for i in range(1,23)]+["chrX"])
    return A[:,ok].copy()

# common genes across adata & bm
common=adata.var_names.intersection(bm.var_names)
adata=adata[:,common].copy(); bm=bm[:,common].copy()

# pooled reference = BM1-4 + all D0 normal lymphocytes
lymph=["T","CTL","NK","B","ProB","Plasma"]
lymph_ref=adata[(adata.obs.CellType.isin(lymph))&(adata.obs.malignancy=="normal")].copy()
lymph_ref.obs["refsample"]="Dlymph"
print("BM ref cells:",bm.n_obs,"lymph ref cells:",lymph_ref.n_obs)

MALIG_MIN=100
patients=[p for p in adata.obs.patient.unique()
          if (adata.obs.patient.eq(p)&adata.obs.malignancy.eq("malignant")).sum()>=MALIG_MIN]
print("patients:",len(patients),patients)

clone_rows=[]; summ_rows=[]; adata.obs["clone"]="NA"
for p in patients:
    pc=adata[(adata.obs.patient==p)].copy()
    pc.obs["role"]="obs"
    r1=bm.copy(); r1.obs["role"]="ref"; r1.obs["patient"]=p
    r2=lymph_ref[lymph_ref.obs.patient!=p].copy(); r2.obs["role"]="ref"  # other-patient lymphs
    for c in ["CellType","malignancy","sample","MutTranscripts","WtTranscripts","CyclingScore"]:
        if c not in r1.obs: r1.obs[c]="ref"
    sub=ad.concat([pc,r1,r2],join="inner")
    sub=prep(sub)
    # normalize
    sc.pp.normalize_total(sub,target_sum=1e4); sc.pp.log1p(sub)
    cnv.tl.infercnv(sub,reference_key="role",reference_cat=["ref"],window_size=100,step=10)
    obs=sub[sub.obs.role=="obs"].copy()
    cnv.tl.pca(obs,use_rep="cnv",key_added="cnv_pca",n_comps=min(30,obs.n_obs-1))
    cnv.pp.neighbors(obs,use_rep="cnv_pca")
    cnv.tl.leiden(obs,resolution=0.4)
    Xc=sub.obsm["X_cnv"]
    mag=np.asarray((Xc.multiply(Xc)).mean(1)).ravel()
    sub.obs["cnv_mag"]=mag
    ref_mag=sub.obs.loc[sub.obs.role=="ref","cnv_mag"].mean()
    obs.obs["cnv_mag"]=sub.obs.loc[obs.obs_names,"cnv_mag"].values
    for g,idx in obs.obs.groupby("cnv_leiden").groups.items():
        d=obs.obs.loc[idx]; n=len(d); fm=(d.malignancy=="malignant").mean()
        ratio=d.cnv_mag.mean()/ref_mag
        isc=(n>=30)and(fm>=0.5)and(ratio>1.5)
        clone_rows.append(dict(patient=p,leiden=g,n=n,frac_malig=round(fm,3),
                               cnv_ratio=round(ratio,3),is_malig_clone=isc))
        if isc:
            cells=d.index[d.malignancy=="malignant"]
            adata.obs.loc[cells,"clone"]=f"{p}_c{g}"
    mm=obs.obs.loc[obs.obs.malignancy=="malignant","cnv_mag"].mean()
    ncl=sum(1 for r in clone_rows if r["patient"]==p and r["is_malig_clone"])
    summ_rows.append(dict(patient=p,n_malig=(obs.obs.malignancy=="malignant").sum(),
                          n_ref=int((sub.obs.role=="ref").sum()),
                          malig_cnv_ratio=round(mm/ref_mag,3),n_clones=ncl))
    print(f"{p}: malig_ratio={mm/ref_mag:.2f} clones={ncl} nclusters={obs.obs.cnv_leiden.nunique()}")

pd.DataFrame(clone_rows).to_csv("/tmp/aml/aml_clusters_all.csv",index=False)
pd.DataFrame(summ_rows).to_csv("/tmp/aml/aml_infercnv_summary.csv",index=False)
# store cnv_pca + clone in adata for validation; recompute obs cnv per patient not kept globally
adata.write("/tmp/aml/aml_cnv.h5ad")
print(pd.DataFrame(summ_rows).to_string())
print("total clone cells:",(adata.obs.clone!="NA").sum(),"n clones:",adata.obs.clone.nunique()-1)
