"""
headtohead_repr.py — the untested-core experiment as a reusable, powered-when-liu-lands test.

Benchmarks three predictors of ICB response under BOTH random 5-fold CV and the honest
leave-one-cohort-out (LOCO) scheme, with a LOCO permutation null:
  - immune_floor : IFN / T-cell-inflamed signature score (computed inside each fold)
  - learned_rep  : unsupervised PCA embedding of top-variable genes (fit inside each fold)
  - floor+rep    : both

At n=40 (gide+riaz) the verdict was: immune floor LOCO AUROC 0.70 (perm p=0.004);
learned rep 0.56 (p=0.16, does NOT beat chance across cohorts); combining DEGRADES the floor.
The random-CV parity (0.75 vs 0.76) was within-cohort overfitting (rep variance 6x the floor).
When liu2019 lands (n~132, 3 cohorts) this re-runs and the LOCO test gains a third held-out
cohort — the decisive power increase for whether a learned RNA representation carries any
cross-cohort ICB signal beyond immune composition.
"""
from __future__ import annotations
import json, urllib.request
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score

IFN=["IFNG","STAT1","CCR5","CXCL9","CXCL10","CXCL11","IDO1","PRF1","GZMA","MHC2TA","HLA-DRA",
     "CD3D","CD3E","CD2","IL2RG","NKG7","CCL5","LAG3","TAGAP","CD8A","CD27","CXCR6"]

def _ens(symbols):
    req=urllib.request.Request("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=json.dumps({"symbols":symbols}).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"})
    r=json.load(urllib.request.urlopen(req,timeout=90))
    return {s:i["id"] for s,i in r.items() if isinstance(i,dict) and str(i.get("id","")).startswith("ENSG")}

def run(logX, y, cohort, k=10, ntop=2000, n_perm=500, seed=0):
    """logX: log2(TPM+1), rows=samples, cols=ENSG. y: 0/1 response. cohort: group labels."""
    y=np.asarray(y).astype(int); cohort=np.asarray(cohort)
    ifn_ensg=[e for e in _ens(IFN).values() if e in logX.columns]
    def floor(tr,te):
        mu,sd=logX.iloc[tr][ifn_ensg].mean(),logX.iloc[tr][ifn_ensg].std().replace(0,1)
        return (((logX.iloc[tr][ifn_ensg]-mu)/sd).mean(axis=1).values.reshape(-1,1),
                ((logX.iloc[te][ifn_ensg]-mu)/sd).mean(axis=1).values.reshape(-1,1))
    def rep(tr,te):
        L=logX.iloc[tr]; expr=(L>1).mean(axis=0)>0.5; Xe=L.loc[:,expr]
        top=Xe.var(axis=0).sort_values(ascending=False).index[:ntop]
        scl=StandardScaler().fit(L[top]); pca=PCA(n_components=min(k,len(tr)-2),random_state=seed).fit(scl.transform(L[top]))
        return pca.transform(scl.transform(logX.iloc[tr][top])), pca.transform(scl.transform(logX.iloc[te][top]))
    def both(tr,te):
        a=floor(tr,te); b=rep(tr,te); return np.hstack([a[0],b[0]]),np.hstack([a[1],b[1]])
    builders={"immune_floor":floor,"learned_rep":rep,"floor+rep":both}
    def auroc(build,splitter,groups=None,yy=None):
        yy=y if yy is None else yy; oof=np.full(len(yy),np.nan); idx=np.arange(len(yy))
        it=splitter.split(idx,yy,groups) if groups is not None else splitter.split(idx,yy)
        for tr,te in it:
            if len(np.unique(yy[tr]))<2: continue
            Xtr,Xte=build(tr,te)
            oof[te]=LogisticRegression(max_iter=2000,C=0.5).fit(Xtr,yy[tr]).predict_proba(Xte)[:,1]
        m=~np.isnan(oof)
        return roc_auc_score(yy[m],oof[m]) if len(np.unique(yy[m]))>1 else np.nan
    skf=StratifiedKFold(5,shuffle=True,random_state=seed); logo=LeaveOneGroupOut()
    ncoh=len(np.unique(cohort))
    out={"n":int(len(y)),"n_cohorts":int(ncoh),"cv5":{},"loco":{},"loco_perm_p":{}}
    for nm,b in builders.items(): out["cv5"][nm]=float(auroc(b,skf))
    if ncoh>=2:
        rng=np.random.default_rng(seed)
        obs={nm:float(auroc(b,logo,groups=cohort)) for nm,b in builders.items()}
        out["loco"]=obs
        null={nm:[] for nm in builders}
        for _ in range(n_perm):
            yp=rng.permutation(y)
            for nm,b in builders.items():
                v=auroc(b,logo,groups=cohort,yy=yp)
                if not np.isnan(v): null[nm].append(v)
        for nm in builders:
            nd=np.array(null[nm]); out["loco_perm_p"][nm]=float((nd>=obs[nm]).mean()) if len(nd) else None
    return out
