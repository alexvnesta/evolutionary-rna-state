"""Evo2-7B scoring of novel-junction sequences (within-Gide encoder pass).
For each junction: mean log-likelihood of spliced vs contiguous reference sequence.
delta = LL(spliced) - LL(contiguous) = how much LESS plausible the aberrant splice is.
Batched, checkpoints to out/ periodically."""
import json, os, time, torch
os.makedirs("out", exist_ok=True)
from evo2 import Evo2
t0=time.time()
m=Evo2("evo2_7b")
print("model loaded", round(time.time()-t0), flush=True)
recs=json.load(open("junction_seqs.json"))
print("scoring", len(recs), "junctions x2 seqs", flush=True)

@torch.inference_mode()
def ll(seq):
    ids=torch.tensor(m.tokenizer.tokenize(seq),dtype=torch.int).unsqueeze(0).to("cuda:0")
    out=m.model.forward(ids)
    logits=out[0] if isinstance(out,(tuple,list)) else out
    logp=torch.log_softmax(logits.float(),dim=-1)
    tgt=ids[:,1:]
    lp=logp[:,:-1,:].gather(-1,tgt.long().unsqueeze(-1)).squeeze(-1)
    return float(lp.mean().item())

res=[]
for i,r in enumerate(recs):
    s=ll(r["spliced"]); c=ll(r["contiguous"])
    res.append({"jid":r["jid"],"chrom":r["chrom"],"istart":r["istart"],"iend":r["iend"],
                "strand":r["strand"],"ll_spliced":s,"ll_contig":c,"delta":s-c})
    if (i+1)%200==0:
        json.dump(res, open("out/evo2_junction_scores.json","w"))
        print(f"  {i+1}/{len(recs)} elapsed {round(time.time()-t0)}s", flush=True)
json.dump(res, open("out/evo2_junction_scores.json","w"))
print("DONE", len(res), "wall", round(time.time()-t0), flush=True)
