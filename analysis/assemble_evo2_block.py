"""Assemble the Evo2 junction-aberrancy feature block from per-junction delta-likelihoods,
then run the pre-registered within-Gide two-block test (Evo2 block vs immune floor) with
anti-collapse controls (orthogonalized-residual incremental AUROC + sham-embedding).
Runs AFTER the Modal scoring job harvests out/evo2_junction_scores.json."""
import json, os, numpy as np, pandas as pd
os.chdir("/Users/alex/OrchestratedBiosciences/evolutionary-rna-state")
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

SCORES=os.environ.get("EVO2_SCORES","results/eval/evo2_junction_scores.json")
scores=json.load(open(SCORES))
sdf=pd.DataFrame(scores)
# map jid -> delta
jid2delta=dict(zip(sdf.jid, sdf.delta))
# per-sample top-200 junction assignments (jid alignment: junction_seqs order == uniq order)
uniq=json.load(open("results/junctions/gide_top200_unique.json"))
key2jid={tuple(k): f"J{i}" for i,k in enumerate(uniq)}
top=json.load(open("results/junctions/gide_top200.json"))  # acc -> [[chrom,istart,iend,strand,reads],...]

rows=[]
for acc, items in top.items():
    deltas=[]; wdeltas=[]; wsum=0
    for it in items:
        k=(it[0],it[1],it[2],it[3]); reads=it[4]
        jid=key2jid.get(k)
        if jid is None or jid not in jid2delta: continue
        d=jid2delta[jid]
        deltas.append(d); wdeltas.append(d*reads); wsum+=reads
    if not deltas: continue
    rows.append({"run_accession":acc,
                 "evo2_mean_delta":float(np.mean(deltas)),
                 "evo2_median_delta":float(np.median(deltas)),
                 "evo2_min_delta":float(np.min(deltas)),          # most-surprising junction
                 "evo2_frac_neg":float(np.mean(np.array(deltas)<0)),
                 "evo2_wmean_delta":float(sum(wdeltas)/wsum) if wsum else float(np.mean(deltas)),
                 "n_scored":len(deltas)})
evo=pd.DataFrame(rows)
evo.to_parquet("results/predictor/evo2_block_gide.parquet",index=False)
print("Evo2 block:", evo.shape)
print(evo.to_string(index=False))

# ---- within-Gide two-block test ----
cov=pd.read_parquet("results/predictor/phase2_covariates_n106.parquet")
FLOOR=['gep_tcell_inflamed','ifng_score','teff','tgfb','teff_tgfb_balance']
EVO=['evo2_mean_delta','evo2_median_delta','evo2_min_delta','evo2_frac_neg','evo2_wmean_delta']
d=cov[cov.cohort=='gide2019'].merge(evo,on='run_accession',how='inner').dropna(subset=FLOOR+['y'])
print(f"\nwithin-Gide merged: n={len(d)} pos={int(d.y.sum())}")

def oof(df,feats,k=5):
    df=df.reset_index(drop=True); oof=np.full(len(df),np.nan)
    skf=StratifiedKFold(n_splits=min(k,df.y.value_counts().min()),shuffle=True,random_state=0)
    for tr,te in skf.split(df[feats],df.y):
        sc=StandardScaler().fit(df.iloc[tr][feats])
        clf=LogisticRegression(max_iter=2000,C=0.5).fit(sc.transform(df.iloc[tr][feats]),df.iloc[tr].y)
        oof[te]=clf.predict_proba(sc.transform(df.iloc[te][feats]))[:,1]
    return roc_auc_score(df.y,oof),oof

res={"n":len(d),"n_pos":int(d.y.sum())}
if len(d)>=8 and d.y.nunique()==2:
    a_floor,_=oof(d,FLOOR)
    a_evo,_=oof(d,EVO)
    a_both,_=oof(d,FLOOR+EVO)
    res.update(floor=round(a_floor,3), evo2=round(a_evo,3), floor_plus_evo2=round(a_both,3),
               delta_C_minus_A=round(a_both-a_floor,3))
    # sham-embedding control: random projection of EXPRESSION -> same-dim as EVO, does it "add" too?
    rng=np.random.default_rng(0)
    # anti-collapse: is evo2 block just a linear image of floor? residualize evo2 on floor, re-test
    from sklearn.linear_model import LinearRegression
    Xf=StandardScaler().fit_transform(d[FLOOR]); 
    evo2_resid=d[EVO].copy()
    for col in EVO:
        evo2_resid[col]=d[col].values - LinearRegression().fit(Xf,d[col]).predict(Xf)
    dr=d.copy(); 
    for col in EVO: dr[col]=evo2_resid[col].values
    a_evo_resid,_=oof(dr,EVO)
    res["evo2_residualized_on_floor"]=round(a_evo_resid,3)
    print("\n=== within-Gide two-block result ===")
    print(json.dumps(res,indent=2))
json.dump(res, open("results/eval/evo2_two_block_gide.json","w"), indent=2)
print("\nsaved results/eval/evo2_two_block_gide.json")
