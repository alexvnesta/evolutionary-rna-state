"""ICB outcome test on Riaz melanoma bulk (GSE91061).
Scores the malignant-intrinsic IFN/MHC signature and an immune-infiltration proxy in
pre-treatment bulk, tests each against RECIST response, and checks whether IFN/MHC adds
beyond infiltration. Inputs: GSE91061 FPKM (Entrez-indexed) + per-sample response labels.
"""
import numpy as np, pandas as pd, json
from scipy.stats import mannwhitneyu, pearsonr
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm

ifn_entrez=[3105,3106,3107,567,6772,3659,6890,6891,2633,3627,5698,3133,84166]  # HLA-A/B/C,B2M,STAT1,IRF1,TAP1/2,GBP1,CXCL10,PSMB9,HLA-E,NLRC5
inf_entrez=[5788,915,916,925,968,4069,931,4818,914,3684,3002,3001]              # PTPRC,CD3D/E,CD8A,CD68,LYZ,MS4A1,NKG7,CD2,ITGAM,GZMB/A

fpkm=pd.read_csv("fpkm.csv.gz",index_col=0); fpkm.index=fpkm.index.astype(int)
L=np.log2(fpkm+1); Z=L.sub(L.mean(1),axis=0).div(L.std(1)+1e-9,axis=0)
sc=lambda e:Z.loc[[g for g in e if g in Z.index]].mean(0)
ifn,inf=sc(ifn_entrez),sc(inf_entrez)
pre=json.load(open("riaz_pre.json"))["pre"]; lab={x["title"]:x["bin"] for x in pre if x["bin"]}
cols=[c for c in fpkm.columns if c in lab]; y=np.array([1 if lab[c]=="R" else 0 for c in cols])
a,b=ifn[cols].values,inf[cols].values
for n,s in [("IFN/MHC",a),("infiltration",b)]:
    print(n,"AUC=%.3f MWU_p=%.3f"%(roc_auc_score(y,s),mannwhitneyu(s[y==1],s[y==0])[1]))
print("IFN vs inf r=%.2f"%pearsonr(a,b)[0])
m=sm.Logit(y,sm.add_constant(np.column_stack([b,a]))).fit(disp=0)
print("joint logit: inf_p=%.3f ifn_p=%.3f"%(m.pvalues[1],m.pvalues[2]))
