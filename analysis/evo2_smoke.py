"""Minimal Evo2-7B GPU smoke test: load from hydrated weights volume, score 2 short sequences.
Proves the encoder path end-to-end. Writes a small JSON to out/."""
import json, os, time, pathlib
os.makedirs("out", exist_ok=True)
t0=time.time()
import torch
from evo2 import Evo2
print("torch", torch.__version__, "| cuda", torch.cuda.is_available(), "| dev", torch.cuda.get_device_name() if torch.cuda.is_available() else "none", flush=True)
m = Evo2("evo2_7b")
print("model loaded in", round(time.time()-t0), "s", flush=True)

# two short DNA sequences: a plausible one and a scrambled one — delta-likelihood sanity
seqs = {
  "ref_like": "ATGGCCTGCACTGGAAGCTTCAAGCTGACCGTGGACGGCACCAAGTTCGAGGTGAAG",
  "scrambled": "TGACGTACGTTGCAAGCTAGCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGC",
}
def score(seq):
    ids = torch.tensor(m.tokenizer.tokenize(seq), dtype=torch.int).unsqueeze(0).to("cuda:0")
    with torch.inference_mode():
        out = m.model.forward(ids)
        logits = out[0] if isinstance(out,(tuple,list)) else out
        logp = torch.log_softmax(logits.float(), dim=-1)
        # mean log-likelihood of the realized tokens (shifted)
        tgt = ids[:,1:]; lp = logp[:,:-1,:].gather(-1, tgt.long().unsqueeze(-1)).squeeze(-1)
        return float(lp.mean().item())
res = {k: score(v) for k,v in seqs.items()}
res["wall_s"] = round(time.time()-t0,1)
res["gpu"] = torch.cuda.get_device_name()
print("SCORES:", json.dumps(res), flush=True)
json.dump(res, open("out/evo2_smoke.json","w"), indent=2)
