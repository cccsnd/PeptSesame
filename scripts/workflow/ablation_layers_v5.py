#!/usr/bin/env python3
"""Evidence-layer contribution analysis (M1-M4) for PeptSesame.

M1 motif only → M2 +conservation → M3 +structure → M4 +expression → M5 full
Metrics: Spearman rho vs full ranking, top-50/top-100 Jaccard overlap.
Requires a completed Layer 2 run (results/02_layer2_scoring/<species>/scored_sorfs.tsv).
Output: results/08_benchmark/ablation_layers.tsv under PEPTSESAME_ROOT.
"""
import csv
import math
import os
from pathlib import Path

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("PEPTSESAME_ROOT", os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
V5 = Path(ROOT)
OUT = V5 / "results/08_benchmark/ablation_layers.tsv"

CH = ["sequence_features", "conservation", "expression", "structural", "motif"]
W = {"sequence_features": 0.20, "conservation": 0.15, "expression": 0.10,
     "structural": 0.15, "motif": 0.15, "ml": 0.25}


def load():
    rows = []
    with open(V5 / "02_layer2_scoring/Yu11/scored_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append({c: float(r[c]) for c in CH})
    return rows


def score(rows, chans, wts):
    return [sum(wts[c] * r[c] for c in chans) for r in rows]


def rho(a, b):
    n = len(a)
    ra = [0] * n
    for r, i in enumerate(sorted(range(n), key=lambda i: a[i], reverse=True)):
        ra[i] = r
    rb = [0] * n
    for r, i in enumerate(sorted(range(n), key=lambda i: b[i], reverse=True)):
        rb[i] = r
    ma = sum(ra) / n
    mb = sum(rb) / n
    cov = sum((ra[i]-ma)*(rb[i]-mb) for i in range(n))
    va = math.sqrt(sum((x-ma)**2 for x in ra))
    vb = math.sqrt(sum((x-mb)**2 for x in rb))
    return cov / (va*vb) if va and vb else 0.0


def jac(a, b, k):
    ta = set(sorted(range(len(a)), key=lambda i: a[i], reverse=True)[:k])
    tb = set(sorted(range(len(b)), key=lambda i: b[i], reverse=True)[:k])
    return len(ta & tb) / len(ta | tb) if ta | tb else 0.0


def main():
    rows = load()
    n = len(rows)
    print(f"Yu11 候选: {n:,}")
    full = score(rows, CH, W)
    models = {
        "M1 motif only": (["motif"], W),
        "M2 + conservation": (["motif", "conservation"], W),
        "M3 + structure": (["motif", "conservation", "structural"], W),
        "M4 + expression": (["motif", "conservation", "structural", "expression"], W),
        "M5 full (5 channels + neutral ML)": (CH, W),
    }
    out = [("model", "spearman_vs_full", "top50_jaccard", "top100_jaccard")]
    for name, (chans, wts) in models.items():
        s = score(rows, chans, wts)
        r = rho(s, full)
        j50 = jac(s, full, 50)
        j100 = jac(s, full, 100)
        print(f"{name:<32} ρ={r:.4f}  top50 J={j50:.3f}  top100 J={j100:.3f}")
        out.append((name, f"{r:.4f}", f"{j50:.3f}", f"{j100:.3f}"))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerows(out)
    print(f"✅ {OUT}")


if __name__ == "__main__":
    main()
