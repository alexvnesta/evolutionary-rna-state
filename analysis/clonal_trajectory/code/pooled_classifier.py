"""
Pooled powered supervised classifier across cSCC + DDLPS + AML.

Target: a HARMONIZED immune-state label (immune-hot vs immune-cold), the shared axis
that recurred in all three cohorts. Raw class_k2 is NOT consistently oriented across
cohorts (cSCC hot=class2; DDLPS hot=class1; AML hot=class1), so we orient each cohort's
label by its interferon_MHC program mean.

Design choices to avoid leakage / cohort-shortcut:
- Features standardized WITHIN cohort (z-score per gene per cohort) so absolute
  tissue-baseline offsets cannot be used as a shortcut.
- Cross-validation:
    (1) Leave-one-cohort-out (LOCO): strongest test of transferability.
    (2) Grouped-by-patient StratifiedGroupKFold: pooled within-cohort generalization.
- Permutation null (labels shuffled within cohort) for significance.
- Connect balanced accuracy back to the calibrated power grid (N vs effect size).
"""
import numpy as np, pandas as pd, os, json, host
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

TAB="/Users/alex/OrchestratedBiosciences/evolutionary-rna-state/analysis/clonal_trajectory/tables"

def load_cohort(name, pb_npy, meta_csv, genes_csv, class_df, ifn_col):
    PB=np.load(pb_npy)
    meta=pd.read_csv(meta_csv)
    genes=pd.read_csv(genes_csv).iloc[:,0].astype(str).values
    # attach class if meta lacks it
    if "class_k2" not in meta.columns:
        meta=meta.merge(class_df[["clone","class_k2"]], on="clone", how="left")
    # orient immune-hot: the class with higher mean IFN/MHC -> label 1 (hot)
    g=class_df.groupby("class_k2")[ifn_col].mean()
    hot=g.idxmax()
    meta["immune_hot"]=(meta["class_k2"]==hot).astype(int)
    meta=meta.dropna(subset=["class_k2"]).reset_index(drop=True)
    PB=PB[meta.index.values] if len(meta)==PB.shape[0] else PB  # aligned
    meta["cohort"]=name
    meta["group"]=name+"::"+meta["patient"].astype(str)
    return PB, meta, genes

# --- load three cohorts ---
cscc_cls=pd.read_csv(host.artifact_path("1fc5f377-4083-427d-a924-e9d69d57b9a7"))
ddlps_cls=pd.read_csv(host.artifact_path("e1773bac-3f0c-4e09-9079-3ebc0ef5380b"))
aml_cls=pd.read_csv(host.artifact_path("2cd3adcb-6012-4dce-aac6-d6f6240018ab"))

PBc,Mc,Gc=load_cohort("cSCC","/tmp/cscc/clone_pseudobulk.npy","/tmp/cscc/clone_pseudobulk_meta.csv",
                      "/tmp/cscc/pb_genes.csv", cscc_cls, "interferon_MHC")
PBd,Md,Gd=load_cohort("DDLPS","/tmp/ddlps/ddlps_clone_pseudobulk.npy","/tmp/ddlps/ddlps_clone_pseudobulk_meta.csv",
                      "/tmp/ddlps/ddlps_pb_genes.csv", ddlps_cls, "interferon_MHC")
PBa,Ma,Ga=load_cohort("AML", f"{TAB}/aml_clone_pseudobulk.npy", f"{TAB}/aml_clone_pseudobulk_meta.csv",
                      f"{TAB}/aml_pb_genes.csv", aml_cls, "interferon_MHC")

# --- shared genes ---
shared=sorted(set(Gc)&set(Gd)&set(Ga))
print("shared genes:",len(shared))
def sub(PB,G): 
    idx={g:i for i,g in enumerate(G)}; return PB[:,[idx[g] for g in shared]]
Xc,Xd,Xa=sub(PBc,Gc),sub(PBd,Gd),sub(PBa,Ga)

# --- within-cohort z-score standardization ---
def zc(X): return StandardScaler().fit_transform(X)
Xc,Xd,Xa=zc(Xc),zc(Xd),zc(Xa)

X=np.vstack([Xc,Xd,Xa])
meta=pd.concat([Mc,Md,Ma],ignore_index=True)
y=meta["immune_hot"].values
groups=meta["group"].values
cohort=meta["cohort"].values
print("pooled X:",X.shape,"| y balance:",dict(pd.Series(y).value_counts()))
print("by cohort:\n",meta.groupby("cohort")["immune_hot"].agg(["count","mean"]).to_string())

# feature selection: top-variance genes (on pooled z-scored) to keep p manageable
K=2000
v=X.var(0); topk=np.argsort(v)[::-1][:K]
Xk=X[:,topk]

def fit_eval(Xtr,ytr,Xte):
    clf=LogisticRegression(penalty="l2",C=0.05,class_weight="balanced",max_iter=2000)
    clf.fit(Xtr,ytr); return clf.predict(Xte), clf.decision_function(Xte)

# --- (1) Leave-one-cohort-out ---
loco={}
for held in ["cSCC","DDLPS","AML"]:
    tr=cohort!=held; te=cohort==held
    pred,score=fit_eval(Xk[tr],y[tr],Xk[te])
    ba=balanced_accuracy_score(y[te],pred)
    try: auc=roc_auc_score(y[te],score)
    except: auc=np.nan
    loco[held]={"n_test":int(te.sum()),"bacc":round(ba,3),"auc":round(auc,3)}
print("\nLOCO:",json.dumps(loco,indent=1))

# --- (2) grouped-by-patient StratifiedGroupKFold ---
sgkf=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=0)
preds=np.zeros_like(y); scores=np.zeros(len(y),float)
for tr,te in sgkf.split(Xk,y,groups):
    p,s=fit_eval(Xk[tr],y[tr],Xk[te]); preds[te]=p; scores[te]=s
ba_cv=balanced_accuracy_score(y,preds); auc_cv=roc_auc_score(y,scores)
print(f"\nGrouped 5-fold: bacc={ba_cv:.3f} auc={auc_cv:.3f} (N={len(y)})")

# --- permutation null (shuffle labels within cohort) ---
rng=np.random.default_rng(0); perm_ba=[]
for _ in range(200):
    yp=y.copy()
    for c in np.unique(cohort):
        m=cohort==c; yp[m]=rng.permutation(yp[m])
    pr=np.zeros_like(yp)
    for tr,te in sgkf.split(Xk,yp,groups):
        pr[te]=fit_eval(Xk[tr],yp[tr],Xk[te])[0]
    perm_ba.append(balanced_accuracy_score(yp,pr))
perm_ba=np.array(perm_ba); pval=(np.sum(perm_ba>=ba_cv)+1)/(len(perm_ba)+1)
print(f"perm null bacc={perm_ba.mean():.3f}+/-{perm_ba.std():.3f}; p={pval:.4f}")

# save results
out={"n_clones":int(len(y)),"n_shared_genes":len(shared),"k_features":K,
     "y_balance":{int(k):int(v) for k,v in pd.Series(y).value_counts().items()},
     "by_cohort":meta.groupby("cohort")["immune_hot"].agg(["count","mean"]).round(3).to_dict(),
     "loco":loco,"grouped_cv_bacc":round(ba_cv,3),"grouped_cv_auc":round(auc_cv,3),
     "perm_null_mean":round(float(perm_ba.mean()),3),"perm_null_sd":round(float(perm_ba.std()),3),
     "perm_p":round(float(pval),4)}
json.dump(out,open("/tmp/pooled_results.json","w"),indent=1)
meta_out=meta[["cohort","patient","clone","class_k2","immune_hot"]].copy()
meta_out["cv_pred"]=preds; meta_out["cv_score"]=scores.round(3)
meta_out.to_csv("/tmp/pooled_labels_preds.csv",index=False)
print("\nSAVED /tmp/pooled_results.json and /tmp/pooled_labels_preds.csv")
print(json.dumps(out,indent=1))
