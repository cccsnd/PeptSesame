#!/usr/bin/env python3
"""benchmark_weight_sensitivity.py — 评分权重敏感性 (08_benchmark, current pipeline)

Method (consistent with the earlier score_weight_sensitivity.py; 5→6 channels):
- 输入: results/02_layer2_scoring/Yu11/scored_sorfs.tsv (445 万 sORF, 513MB)
- 权重 (pipeline/layer2_scoring/scoring_core.py DEFAULT_WEIGHTS, 和≈1.0):
    sequence_features 0.20, ml 0.25, conservation 0.15,
    expression 0.10, structural 0.15, motif 0.15
- 扰动: 每通道 ±20% / ±50% (其余通道权重不变, 归一化后重算 aggregated)
- 指标: ① 全量 aggregated_score Spearman 秩相关 vs 原始
        ② Top-100 候选 Jaccard 重叠率
- 结论判断: ρ>0.98 且 Top-100 Jaccard 大多 >0.7 → 排名对权重稳健

输出: results/08_benchmark/weight_sensitivity.tsv + 控制台摘要
"""
import csv
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results"
SRC = R2 / "02_layer2_scoring/Yu11/scored_sorfs.tsv"
OUT_DIR = R2 / "08_benchmark"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WEIGHTS = {
    "sequence_features": 0.20,
    "ml": 0.25,
    "conservation": 0.15,
    "expression": 0.10,
    "structural": 0.15,
    "motif": 0.15,
}
# scored_sorfs.tsv 列名 -> 内部权重键
COL_MAP = {
    "sequence_features": "sequence_features",
    "conservation": "conservation",
    "expression": "expression",
    "structural": "structural",
    "motif": "motif",
    "ml_score": "ml",
}
PERTURB = [0.2, 0.5]


def main():
    print(f"读取 {SRC} ...")
    df = pd.read_csv(SRC, sep="\t", usecols=list(COL_MAP.keys()))
    print(f"  {len(df):,} sORF × {len(COL_MAP)} 通道")

    X = df[list(COL_MAP.keys())].to_numpy(dtype=np.float64)
    w0 = np.array([WEIGHTS[COL_MAP[c]] for c in COL_MAP])
    w0 = w0 / w0.sum()

    # 原始 aggregated
    agg0 = X @ w0

    # Top-100 (原始)
    top0 = set(np.argsort(-agg0)[:100].tolist())

    rows = []
    for col_name, key in COL_MAP.items():
        j = list(COL_MAP.keys()).index(col_name)
        for frac in PERTURB:
            for sign, tag in [(1, "+"), (-1, "-")]:
                w = w0.copy()
                w[j] *= (1 + sign * frac)
                w = w / w.sum()
                agg = X @ w
                rho = spearmanr(agg0, agg).statistic
                top1 = set(np.argsort(-agg)[:100].tolist())
                jac = len(top0 & top1) / len(top0 | top1)
                rows.append({
                    "channel": key, "perturbation": f"{tag}{frac:.0%}",
                    "spearman_rho": round(float(rho), 4),
                    "top100_jaccard": round(jac, 3),
                })
                print(f"  {key:16s} {tag}{frac:.0%}: ρ={rho:.4f}  Jaccard={jac:.3f}")

    # 输出
    with open(OUT_DIR / "weight_sensitivity.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["channel", "perturbation",
                                       "spearman_rho", "top100_jaccard"])
        w.writeheader()
        w.writerows(rows)

    rho_min = min(r["spearman_rho"] for r in rows)
    jac_min = min(r["top100_jaccard"] for r in rows)
    print(f"\n✅ weight_sensitivity.tsv  (ρ min={rho_min}, Jaccard min={jac_min})")
    print("结论: ρ 全 > 0.98 = 排名稳健; 最敏感通道看 Jaccard 最低者 (通常 ml)")


if __name__ == "__main__":
    main()
