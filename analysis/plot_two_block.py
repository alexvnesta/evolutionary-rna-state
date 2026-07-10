#!/usr/bin/env python3
"""Two-block result figure: floor vs non-ref vs combined AUROC, with CI + perm-p + delta.
Consumes results/eval/nonref_vs_floor_*.json. Publication-grade per figure-style skill.
"""
import json, sys, os
import numpy as np, matplotlib.pyplot as plt

ERS = "/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"

def make_figure(result_json, out_png):
    R = json.load(open(result_json))
    blocks = R["blocks"]
    order = [("A_floor", "Immune floor\n(IFN-\u03b3 / T-cell, 5 feat)"),
             ("B_nonref", f"Non-reference\n(editing/IR/splice/TE, {R['nonref_n_features']} feat)"),
             ("C_floor_plus_nonref", "Floor + non-reference")]
    focal = "#C1272D"; floor_c = "#3B6BA5"; comb_c = "#6B4C9A"
    colors = {"A_floor": floor_c, "B_nonref": focal, "C_floor_plus_nonref": comb_c}
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ys = np.arange(len(order))[::-1]
    for y, (k, lab) in zip(ys, order):
        b = blocks[k]; auc = b["auroc"]
        if auc is None: continue
        lo, hi = b["ci95"]
        ax.plot([lo, hi], [y, y], color=colors[k], lw=2.4, solid_capstyle="round", zorder=2)
        ax.scatter([auc], [y], s=90, color=colors[k], zorder=3, edgecolor="white", linewidth=1.2)
        p = b["perm_p"]
        ptxt = f"perm p={p:.3f}" if p is not None else ""
        ax.text(hi + 0.012, y, f"AUROC={auc:.3f}  {ptxt}", va="center", ha="left", fontsize=7.5)
    ax.axvline(0.5, color="#888888", ls="--", lw=1, zorder=1)
    ax.text(0.5, len(order)-0.35, "chance", color="#888888", fontsize=7, ha="center")
    ax.set_yticks(ys); ax.set_yticklabels([lab for _, lab in order], fontsize=8)
    ax.set_xlabel("Out-of-fold AUROC ("+R["frame"]+")", fontsize=9)
    ax.set_xlim(0.30, 1.02); ax.margins(y=0.15)
    d = R.get("delta_C_minus_A")
    n = R["n"]; cohorts = ", ".join(f"{k}:{v}" for k, v in R["cohorts"].items())
    title = f"Does the non-reference RNA layer add ICB-response signal beyond the immune floor?"
    ax.set_title(title, fontsize=9.5, loc="left", pad=10)
    sub = f"n={n} ({R['n_pos']}R/{R['n_neg']}N) \u00b7 {cohorts}"
    if d is not None: sub += f" \u00b7 \u0394(C\u2212A)={d:+.3f} AUROC"
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=7.5, color="#444444")
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    return R

if __name__ == "__main__":
    j = sys.argv[1] if len(sys.argv) > 1 else f"{ERS}/results/eval/nonref_vs_floor_loco.json"
    o = sys.argv[2] if len(sys.argv) > 2 else f"{ERS}/results/eval/two_block_result.png"
    R = make_figure(j, o); print("wrote", o); print(json.dumps(R, indent=2))
