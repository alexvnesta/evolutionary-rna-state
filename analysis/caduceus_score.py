#!/usr/bin/env python
"""Frozen Caduceus embedding-delta on the SAME junction set the Evo2 encoder test used.
Per junction: mean-pooled last-hidden embedding of the spliced window minus the contiguous
window -> L2 norm of the delta = "how much the splice perturbs the sequence representation".
Aggregated per-sample (over that sample's top-200 junctions) exactly as the Evo2 block was,
so the two-block test is directly comparable to Evo2/HyenaDNA/EVA.

Runs on GPU (mamba-ssm requires CUDA). Reads junction_seqs.json; writes caduceus_delta.json.
"""
import json, sys, time
import numpy as np, torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

MID = "kuleshov-group/caduceus-ps_seqlen-1k_d_model-256_n_layer-4_lr-8e-3"

def main(seqs_path, out_path):
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(MID, trust_remote_code=True)
    m = AutoModelForMaskedLM.from_pretrained(MID, trust_remote_code=True).to(dev).eval()
    seqs = json.load(open(seqs_path))
    print(f"[caduceus] {len(seqs)} junctions on {dev}", flush=True)

    def embed(s):
        ids = tok(s, return_tensors="pt", truncation=True, max_length=1024).to(dev)
        with torch.no_grad():
            out = m(**ids, output_hidden_states=True)
        h = out.hidden_states[-1][0]            # (L, d)
        return h.mean(0).float().cpu().numpy()  # mean-pool -> (d,)

    rows = []
    t0 = time.time()
    for i, j in enumerate(seqs):
        es = embed(j["spliced"]); ec = embed(j["contiguous"])
        d = es - ec
        rows.append({"jid": j["jid"],
                     "coord": f"{j['chrom']}:{j['istart']}-{j['iend']}:{j['strand']}",
                     "delta_l2": float(np.linalg.norm(d)),
                     "delta_cos": float(np.dot(es, ec) / (np.linalg.norm(es)*np.linalg.norm(ec) + 1e-9))})
        if i % 200 == 0:
            print(f"[caduceus] {i}/{len(seqs)} {time.time()-t0:.0f}s", flush=True)
    json.dump(rows, open(out_path, "w"))
    print(f"[caduceus] DONE {len(rows)} in {time.time()-t0:.0f}s -> {out_path}", flush=True)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
