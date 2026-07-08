import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, infercnvpy as cnv
import warnings; warnings.filterwarnings("ignore")
np.random.seed(0)
import infercnvpy.tl._infercnv as _icnv
def _serial_process_map(fn, *iterables, **kw):
    return [fn(*args) for args in zip(*iterables)]
_icnv.process_map = _serial_process_map

a = ad.read_h5ad("/tmp/cscc/cscc_raw.h5ad")

# --- gene positions ---
gp = pd.read_csv("/tmp/cscc/gene_pos.csv")
gp.columns = ["symbol","chr","start","end"]
gp = gp.dropna(subset=["symbol","chr","start"])
# keep canonical chromosomes
keep = [str(i) for i in range(1,23)] + ["X","Y"]
gp = gp[gp["chr"].astype(str).isin(keep)].drop_duplicates("symbol")
gp["chromosome"] = "chr"+gp["chr"].astype(str)
gp = gp.set_index("symbol")

# annotate var
a.var["chromosome"] = gp["chromosome"].reindex(a.var_names).values
a.var["start"] = gp["start"].reindex(a.var_names).values
a.var["end"] = gp["end"].reindex(a.var_names).values
mapped = a.var["chromosome"].notna().sum()
print(f"genes mapped to positions: {mapped}/{a.n_vars}")

# --- normalize (log1p CPM-ish) ---
a.X = a.layers["counts"].copy()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

# reference categories = non-epithelial normal compartments
a.obs["is_epithelial"] = (a.obs["level1_celltype"]=="Epithelial").astype(str)
ref_types = ["Fibroblast","Endothelial Cell","Tcell","Mac","B Cell","NK","CD1C","LC","MDSC","CLEC9A","PDC","ASDC"]
a.obs["cnv_ref"] = np.where(a.obs["level1_celltype"].isin(ref_types), a.obs["level1_celltype"].astype(str), a.obs["level1_celltype"].astype(str))

results = {}
tumor_clone_assign = {}
for pat in sorted(a.obs["patient"].unique()):
    ap = a[a.obs["patient"]==pat].copy()
    # need reference cells present
    refs = [t for t in ref_types if (ap.obs["level1_celltype"]==t).sum()>=20]
    n_epi = int((ap.obs["level1_celltype"]=="Epithelial").sum())
    if not refs or n_epi < 50:
        print(f"{pat}: SKIP (epi={n_epi}, refs={refs})"); continue
    cnv.tl.infercnv(ap, reference_key="level1_celltype", reference_cat=refs, window_size=100, n_jobs=1)
    # PCA + leiden on CNV profile, restricted to epithelial (malignant candidate) cells
    cnv.tl.pca(ap, n_comps=min(20, ap.n_obs-1))
    cnv.pp.neighbors(ap, use_rep="cnv_pca", n_neighbors=15)
    cnv.tl.leiden(ap, resolution=0.4, key_added="cnv_clone")
    # CNV signal per cell = mean abs of smoothed CNV
    cnv_score = np.asarray(np.abs(ap.obsm["X_cnv"].toarray() if hasattr(ap.obsm["X_cnv"],"toarray") else ap.obsm["X_cnv"]).mean(1)).ravel()
    ap.obs["cnv_score"] = cnv_score
    epi = ap.obs["level1_celltype"]=="Epithelial"
    # malignant clones = leiden clusters dominated by epithelial cells with elevated CNV
    clone_tab = pd.crosstab(ap.obs["cnv_clone"], epi)
    epi_score = ap.obs.groupby("cnv_clone")["cnv_score"].mean()
    ref_cnv = ap.obs.loc[~epi,"cnv_score"].mean()
    malig_clones = [c for c in ap.obs["cnv_clone"].cat.categories
                    if clone_tab.loc[c].get(True,0) > clone_tab.loc[c].sum()*0.5
                    and epi_score[c] > 1.5*ref_cnv]
    n_malig_clones = len(malig_clones)
    n_malig_cells = int(ap.obs["cnv_clone"].isin(malig_clones).sum())
    results[pat] = dict(n_cells=ap.n_obs, n_epi=n_epi, n_refs=sum((ap.obs.level1_celltype==t).sum() for t in refs),
                        n_clones_total=ap.obs["cnv_clone"].nunique(), n_malig_clones=n_malig_clones,
                        n_malig_cells=n_malig_cells, ref_cnv=float(ref_cnv), mean_epi_cnv=float(ap.obs.loc[epi,"cnv_score"].mean()))
    # store per-cell clone assignment for malignant cells
    sub = ap.obs.loc[ap.obs["cnv_clone"].isin(malig_clones), ["patient","cnv_clone","cnv_score","level2_celltype"]].copy()
    sub["clone_id"] = pat+"_"+sub["cnv_clone"].astype(str)
    tumor_clone_assign[pat] = sub
    print(f"{pat}: cells={ap.n_obs} epi={n_epi} clones_total={results[pat]['n_clones_total']} "
          f"malig_clones={n_malig_clones} malig_cells={n_malig_cells} "
          f"epiCNV={results[pat]['mean_epi_cnv']:.4f} refCNV={ref_cnv:.4f} ratio={results[pat]['mean_epi_cnv']/ref_cnv:.2f}")

pd.DataFrame(results).T.to_csv("/tmp/cscc/infercnv_summary.csv")
allclones = pd.concat(tumor_clone_assign.values()) if tumor_clone_assign else pd.DataFrame()
allclones.to_csv("/tmp/cscc/tumor_clone_assignments.csv")
print("\nSaved infercnv_summary.csv and tumor_clone_assignments.csv")
print("total malignant cells with clone labels:", len(allclones))
