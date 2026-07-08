import numpy as np, pandas as pd
from scipy.stats import pearsonr, mannwhitneyu
res=pd.read_csv("/tmp/aml/aml_clone_programs_classified.csv")
ifn=res["interferon_MHC"].values
# differentiation axis = myeloid_diff - HSC_progenitor
diffaxis=res["myeloid_diff"].values - res["HSC_progenitor"].values
r,p=pearsonr(diffaxis,ifn)
print(f"differentiation axis (myeloid-HSC) vs IFN/MHC: r={r:.3f} p={p:.4f}")
# class-level difference
c1=ifn[res.class_k2==1]; c2=ifn[res.class_k2==2]
u,pmw=mannwhitneyu(c1,c2,alternative="two-sided")
print(f"IFN/MHC class1({len(c1)}) mean={c1.mean():.3f} vs class2({len(c2)}) mean={c2.mean():.3f}  MWU p={pmw:.4f}")
# save full coupling table
rows=[]
for k in ["HSC_progenitor","myeloid_diff","cell_cycle"]:
    rr,pp=pearsonr(res[k].values,ifn); rows.append(dict(program=k,r_vs_IFN_MHC=round(rr,3),p=round(pp,4)))
rows.append(dict(program="differentiation_axis(myeloid-HSC)",r_vs_IFN_MHC=round(r,3),p=round(p,4)))
rows.append(dict(program="class_IFN_MHC_MWU",r_vs_IFN_MHC=round(c1.mean()-c2.mean(),3),p=round(pmw,4)))
pd.DataFrame(rows).to_csv("/tmp/aml/aml_immune_coupling.csv",index=False)
print(pd.DataFrame(rows).to_string())
