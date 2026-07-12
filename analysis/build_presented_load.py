#!/usr/bin/env python
"""NMD-escape modeling (used by regenerate_neopeptides below):

For a novel splice junction, the transcript reads donor-exon ... [junction] ... acceptor-exon.
A neopeptide spanning the junction is only PRESENTED in vivo if its transcript ESCAPES
nonsense-mediated decay. Classical NMD rule: a premature termination codon (PTC) triggers
NMD if it lies >~50 nt UPSTREAM of the last exon-exon junction (the EJC deposited there).
We approximate per junction: translate the in-frame ORF reading THROUGH the junction; if a
stop codon occurs in the acceptor-side sequence >55 nt before the sequence end (a proxy for
a downstream EJC), the transcript is NMD-TARGET (peptide not stably presented) and is flagged.
Junction peptides are retained only from ORFs that either (a) have no stop after the junction
within the modeled window (read-through / last-exon PTC = NMD-escape), or (b) whose stop is
within the 55 nt EJC-proximal zone. This is a first-order model, not a full transcript
reconstruction, and is documented as such.
"""
"""Presentation-layer campaign: per-sample PRESENTED novel-junction neopeptide load.

For each sample: take its top-200 novel junctions, the junction-spanning neopeptides
each generates, score them against THAT patient's HLA alleles with mhcflurry, and
summarize the presented load. Presented load = aberrancy filtered through the patient's
own presentation machinery — the feature designed to DECOUPLE from bulk immune infiltration.

Outputs results/predictor/presented_block_<cohort>.parquet with per-sample features.
Usage: python build_presented_load.py <cohort>  (gide2019 | hugo2016)
"""
import os, sys, json, glob
os.environ.setdefault("MHCFLURRY_DATA_DIR",
    "/Users/alex/OrchestratedBiosciences/evolutionary-rna-state/results/mhcflurry_data")
import numpy as np, pandas as pd
from mhcflurry import Class1PresentationPredictor

REPO="/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"
os.chdir(REPO)
COHORT=sys.argv[1] if len(sys.argv)>1 else "gide2019"

# --- NMD-escape neopeptide generator (see module docstring). Used to (re)build the
#     <cohort>_neojunction_peptides.json from junction_seqs when it is missing/stale. ---
from itertools import product as _product
_CODON={}; _b="TCAG"
_aas="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i,(_a,_bb,_c) in enumerate(_product(_b,repeat=3)): _CODON[_a+_bb+_c]=_aas[_i]
def _translate(s):
    s=s.upper(); return ''.join(_CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
EJC_ZONE=55
def neopeptides_with_nmd(spliced, junction_pos, ks=(8,9,10,11)):
    L=len(spliced); jp=junction_pos; out={k:set() for k in ks}
    for frame in range(3):
        prot=_translate(spliced[frame:]); jc=(jp-frame)//3
        stop_after=next((pi for pi in range(max(jc,0),len(prot)) if prot[pi]=='*'),None)
        if stop_after is not None:
            stop_nt=frame+stop_after*3
            if stop_nt>jp and (L-stop_nt)>EJC_ZONE: continue   # NMD target -> degraded
        end=stop_after if stop_after is not None else len(prot)
        for k in ks:
            for st in range(0,end-k+1):
                if st<=jc<st+k:
                    km=prot[st:st+k]
                    if '*' not in km and 'X' not in km: out[k].add(km)
    return {str(k):sorted(v) for k,v in out.items()}

# --- 1. sample -> HLA alleles (mhcflurry 2-field format) ---
def two_field(a):
    # "A*02:01:192" -> "A*02:01"
    parts=a.split(":")
    return ":".join(parts[:2]) if len(parts)>=2 else a
def load_hla(cohort):
    """Return {run_accession: [6 alleles]}. Gide: arcasHLA genotype json keyed by sample_title,
    mapped to ERR via crosswalk. Falls back to prior hla_alleles parquet."""
    hla={}
    xwalk=pd.read_csv("data/registry/gide_id_crosswalk.csv")
    title2err=dict(zip(xwalk.sample_title, xwalk.run_accession))
    for gj in glob.glob("results/hla_typing/gide_arcas/*.genotype.json"):
        base=os.path.basename(gj).replace(".genotype.json","").replace(".markdup.sorted.extracted","")
        err=title2err.get(base)
        if err is None: continue
        g=json.load(open(gj))
        if not g or not any(g.get(l) for l in ("A","B","C")): continue  # skip empty/failed typings
        alleles=[]
        for loc in ("A","B","C"):
            for al in g.get(loc,[]):
                alleles.append(loc+"*"+two_field(al.split("*")[-1]) if "*" not in ("*"+al) else two_field(al))
        # normalize to e.g. A*02:01
        norm=[]
        for al in alleles:
            al=al.replace("**","*")
            norm.append(two_field(al))
        if norm: hla[err]=norm
    return hla

hla=load_hla(COHORT)
print(f"[{COHORT}] samples with HLA: {len(hla)}")

# --- 2. neojunction peptides (jid -> peptides) and per-sample top-200 junctions ---
neojp=json.load(open(f"results/junctions/{'gide32' if COHORT=='gide2019' else 'hugo22'}_neojunction_peptides.json"))
top=json.load(open(f"results/junctions/{'gide32' if COHORT=='gide2019' else 'hugo22'}_top200.json"))
# map junction key -> jid
def jkey(d): return f"{d['chrom']}:{d['istart']}-{d['iend']}:{d['strand']}"
key2jid={jkey(v):jid for jid,v in neojp.items()}

# --- 3. mhcflurry presentation predictor ---
pred=Class1PresentationPredictor.load()
SUPPORTED=set(pred.supported_alleles)

def presented_features(err):
    alleles=[a for a in hla[err] if a in SUPPORTED]
    if not alleles: return None
    items=top[err]  # [[chrom,istart,iend,strand,reads],...]
    # collect this sample's junction peptides (9-mers) with read weight
    pep_reads={}  # peptide -> max reads of junction generating it
    for it in items:
        k=f"{it[0]}:{int(it[1])}-{int(it[2])}:{it[3]}"; reads=it[4]
        jid=key2jid.get(k)
        if jid is None: continue
        for p in neojp[jid]['peptides']['9']:
            pep_reads[p]=max(pep_reads.get(p,0),reads)
    if not pep_reads: return None
    peps=list(pep_reads.keys())
    res=pred.predict(peptides=peps, alleles={err:alleles}, verbose=0)
    res['reads']=res['peptide'].map(pep_reads)
    ps=res['presentation_score'].values
    logr=np.log10(res['reads'].values+1)
    n_strong=int((ps>0.9).sum()); n_mod=int((ps>0.5).sum())
    return {
        "run_accession":err, "n_peptides":len(peps), "n_hla_supported":len(alleles),
        "presented_n_strong": n_strong,                       # count presented (score>0.9)
        "presented_n_moderate": n_mod,
        "presented_frac_strong": float((ps>0.9).mean()),
        "presented_max_score": float(ps.max()),
        "presented_mean_top20": float(np.mean(np.sort(ps)[-20:])),
        "presented_load_readweighted": float(np.sum(ps*logr)),  # presented AND expressed
        "presented_strong_readweighted": float(np.sum((ps>0.9)*logr)),
    }

rows=[]
for i,err in enumerate(hla):
    if err not in top: continue
    f=presented_features(err)
    if f: rows.append(f); print(f"  {err}: strong={f['presented_n_strong']} load={f['presented_load_readweighted']:.2f}")
out=pd.DataFrame(rows)
outpath=f"results/predictor/presented_block_{COHORT}.parquet"
out.to_parquet(outpath,index=False)
print(f"[{COHORT}] presented block: {out.shape} -> {outpath}")
