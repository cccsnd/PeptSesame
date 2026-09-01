import os
#!/usr/bin/env python3
"""A. Layer 2 通道消融 (ablation) + 基线比较 — R3/R5 要求

对 Yu11 scored_sorfs.tsv 的五个活跃通道, 比较不同评分方案的排名:
- 单通道: sequence / conservation / expression / structural / motif
- equal-weight (五通道等权, 无 ML)
- 全通道 (生产权重 + ML 中性)
- 基线: motif-only / expression-only / ORF length / AA composition

指标: top-100 Jaccard 重叠 + 全排序 Spearman rho (vs 生产全通道排名)
输出: results/08_benchmark/ablation_v20260825.tsv
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

V5 = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "results"
OUT = V5 / "08_benchmark/ablation_v20260825.tsv"

CHANNELS = ["sequence_features", "conservation", "expression", "structural", "motif"]
WEIGHTS = {"sequence_features": 0.20, "conservation": 0.15, "expression": 0.10,
           "structural": 0.15, "motif": 0.15, "ml": 0.25}
ML_NEUTRAL = 0.5


def load():
    rows = []
    with open(V5 / "02_layer2_scoring/Yu11/scored_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append({c: float(r[c]) for c in CHANNELS})
    return rows


def score_full(r):
    return sum(WEIGHTS[c] * r[c] for c in CHANNELS) + WEIGHTS["ml"] * ML_NEUTRAL


def score_equal(r):
    return sum(r[c] for c in CHANNELS) / len(CHANNELS)


def rank_spearman(a, b):
    """两个分数列表的 Spearman rho (同序数)"""
    n = len(a)
    order_a = sorted(range(n), key=lambda i: a[i], reverse=True)
    rank_a = [0] * n
    for r, i in enumerate(order_a):
        rank_a[i] = r
    order_b = sorted(range(n), key=lambda i: b[i], reverse=True)
    rank_b = [0] * n
    for r, i in enumerate(order_b):
        rank_b[i] = r
    ma = sum(rank_a) / n
    mb = sum(rank_b) / n
    cov = sum((rank_a[i]-ma)*(rank_b[i]-mb) for i in range(n))
    va = math.sqrt(sum((x-ma)**2 for x in rank_a))
    vb = math.sqrt(sum((x-mb)**2 for x in rank_b))
    return cov / (va*vb) if va and vb else 0.0


def top100_jaccard(a, b):
    ta = set(sorted(range(len(a)), key=lambda i: a[i], reverse=True)[:100])
    tb = set(sorted(range(len(b)), key=lambda i: b[i], reverse=True)[:100])
    return len(ta & tb) / len(ta | tb) if ta | tb else 0.0


def main():
    rows = load()
    n = len(rows)
    print(f"Yu11 候选: {n:,}")
    full = [score_full(r) for r in rows]
    schemes = {"full (production weights)": full}
    for c in CHANNELS:
        schemes[f"single:{c}"] = [r[c] for r in rows]
    schemes["equal-weight (5 channels)"] = [score_equal(r) for r in rows]
    # 基线: 长度 (从 seq_id 提取? 用 scored 无长度列 — 用 structural 代理? 改为: 长度不可得, 用 sequence_features 作序列组成基线)
    # scored_sorfs 有 length_aa 列
    with open(V5 / "02_layer2_scoring/Yu11/scored_sorfs.tsv") as f:
        lengths = []
        for r in csv.DictReader(f, delimiter="\t"):
            lengths.append(float(r["length_aa"]))
    schemes["baseline:ORF length"] = lengths

    print(f"{'scheme':<32}{'Spearman vs full':>18}{'Top100 Jaccard':>16}")
    out = [("scheme", "spearman_vs_full", "top100_jaccard_vs_full")]
    for name, s in schemes.items():
        rho = rank_spearman(s, full)
        jac = top100_jaccard(s, full)
        print(f"{name:<32}{rho:>18.4f}{jac:>16.3f}")
        out.append((name, f"{rho:.4f}", f"{jac:.3f}"))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerows(out)
    print(f"✅ {OUT}")


if __name__ == "__main__":
    main()
