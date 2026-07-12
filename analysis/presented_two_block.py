#!/usr/bin/env python
"""Two-block test: does PRESENTED novel-junction neopeptide load add ICB-response
signal beyond the immune floor? Same rigorous protocol as the encoder tests
(20-seed 5-fold CV, fold-contained residualization, cohort-internal permutation).

Also tests presented-load residualized on the RICH immune basis (5 sig + 11 deconv
fractions) — the disentangling frame.
Usage: python presented_two_block.py <cohort>
"""
import os, sys, json
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
REPO="/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"; os.chdir(REPO)
COHORT=sys.argv[1] if len(sys.argv)>1 else "gide2019"

cov=pd.read_parquet("results/predictor/phase2_covariates_n106.parquet")
FLOOR=['gep_tcell_inflamed','ifng_score','teff','tgfb','teff_tgfb_balance']
PRES=['presented_n_strong','presented_frac_strong','presented_max_score',
      'presented_mean_top20','presented_load_readweighted','presented_strong_readweighted']
pres=pd.read_parquet(f"results/predictor/presented_block_{COHORT}.parquet")
try:
    frac=pd.read_parquet(host.artifact_path("8b44598f-94ec-45a9-b702-bd900536478b")) if False else None
except Exception: frac=None
# load deconv fractions from artifact path staged locally
FRACP="results/eval/instaprism_fractions_n106.parquet"
if os.path.exists(FRACP):
    frac=pd.read_parquet(FRACP)
    if 'run_accession' not in frac.columns:
        frac=frac.reset_index().rename(columns={'index':'run_accession'})
    frac=frac.loc[:,~frac.columns.duplicated()]  # drop any duplicated run_accession col
else:
    frac=None
FRAC=[c for c in (frac.columns if frac is not None else []) if c!='run_accession']

d=cov[cov.cohort==COHORT].merge(pres,on='run_accession').dropna(subset=FLOOR+['y']).reset_index(drop=True)
y=d.y.values
def oof_auc(X,y,seed=0):
    oof=np.full(len(y),np.nan); 
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=seed).split(X,y):
        sc=StandardScaler().fit(X[tr]);clf=LogisticRegression(max_iter=2000,C=0.5).fit(sc.transform(X[tr]),y[tr])
        oof[te]=clf.predict_proba(sc.transform(X[te]))[:,1]
    return roc_auc_score(y,oof)
def oof_resid(df,feats,basis,y,seed=0):
    oof=np.full(len(y),np.nan);Xe=df[feats].values;Xf=df[basis].values
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=seed).split(Xe,y):
        scf=StandardScaler().fit(Xf[tr]);Ftr=scf.transform(Xf[tr]);Fte=scf.transform(Xf[te])
        rtr=np.zeros_like(Xe[tr],dtype=float);rte=np.zeros_like(Xe[te],dtype=float)
        for j in range(Xe.shape[1]):
            lr=LinearRegression().fit(Ftr,Xe[tr,j]);rtr[:,j]=Xe[tr,j]-lr.predict(Ftr);rte[:,j]=Xe[te,j]-lr.predict(Fte)
        sc=StandardScaler().fit(rtr);clf=LogisticRegression(max_iter=2000,C=0.5).fit(sc.transform(rtr),y[tr])
        oof[te]=clf.predict_proba(sc.transform(rte))[:,1]
    return roc_auc_score(y,oof)
S=range(20)
res={"cohort":COHORT,"n":int(len(d)),"n_pos":int(y.sum()),
     "floor":round(float(np.mean([oof_auc(d[FLOOR].values,y,s) for s in S])),3),
     "presented_alone":round(float(np.mean([oof_auc(d[PRES].values,y,s) for s in S])),3),
     "floor_plus_presented":round(float(np.mean([oof_auc(d[FLOOR+PRES].values,y,s) for s in S])),3),
     "presented_resid_floor":round(float(np.mean([oof_resid(d,PRES,FLOOR,y,s) for s in S])),3)}
# permutation on the floor-residual
rng=np.random.default_rng(0)
obs=float(np.mean([oof_resid(d,PRES,FLOOR,y,s) for s in range(10)]))
null=np.array([float(np.mean([oof_resid(d,PRES,FLOOR,rng.permutation(y),s) for s in range(3)])) for _ in range(200)])
res["presented_resid_obs10"]=round(obs,3); res["presented_resid_perm_p"]=round(float((np.sum(null>=obs)+1)/201),3)
res["presented_resid_null_mean"]=round(float(null.mean()),3)
np.save(f"results/eval/presented_perm_null_{COHORT}.npy", null)   # real null draws for the figure
# rich-basis disentangling frame (if fractions available)
if frac is not None:
    dr=d.merge(frac,on='run_accession',how='left')
    if dr[FRAC].notna().all(axis=1).sum()==len(dr):
        RICH=FLOOR+FRAC
        res["presented_resid_rich"]=round(float(np.mean([oof_resid(dr,PRES,RICH,y,s) for s in S])),3)
json.dump(res,open(f"results/eval/presented_two_block_{COHORT}.json","w"),indent=2)
print(json.dumps(res,indent=2))
