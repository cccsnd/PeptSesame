import os
#!/usr/bin/env python3
"""去 ML 通道重聚合: 5 通道归一化 (R1/评估要求, 修复稿件-生产不一致)

原 aggregate = 0.20*seq + 0.25*ml + 0.15*cons + 0.10*expr + 0.15*struct + 0.15*motif
新 aggregate = 5 通道归一化: seq/cons/expr/struct/motif 权重 /0.75
  = 0.2667*seq + 0.20*cons + 0.1333*expr + 0.20*struct + 0.20*motif  (范围 0–1)
tier: ≥0.7 high / ≥0.4 medium / <0.4 low (在 0–1 尺度, 不再平移)
对 18 物种 scored_sorfs.tsv 就地重算 aggregated_score + confidence。
"""
import csv
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
V5 = ROOT / "results"
SPECIES = ["Yu11", "14G01", "14G02", "K16", "ken1", "ken8", "S3651",
           "Arabidopsis", "Rice", "Maize", "Soybean", "Tomato", "Grape",
           "Ricinus", "Flax", "Sunflower", "Zisu", "Gastrodia"]

# 5 通道归一化权重 (原 /0.75)
W = {"sequence_features": 0.20/0.75, "conservation": 0.15/0.75,
     "expression": 0.10/0.75, "structural": 0.15/0.75, "motif": 0.15/0.75}
assert abs(sum(W.values()) - 1.0) < 1e-9, sum(W.values())


def tier_of(s):
    return "high" if s >= 0.7 else ("medium" if s >= 0.4 else "low")


summary = []
for sp in SPECIES:
    p = V5 / f"02_layer2_scoring/{sp}/scored_sorfs.tsv"
    if not p.exists():
        print(f"  [skip] {sp}: 无文件")
        continue
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    tiers = {"high": 0, "medium": 0, "low": 0}
    for r in rows:
        s = sum(W[k] * float(r[k]) for k in W)
        r["aggregated_score"] = f"{s:.6f}"
        r["confidence"] = tier_of(s)
        tiers[tier_of(s)] += 1
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    summary.append((sp, n, tiers["high"], tiers["medium"], tiers["low"],
                    f"{tiers['medium']/n*100:.1f}%", f"{tiers['low']/n*100:.1f}%"))
    print(f"  {sp}: n={n:,}  high={tiers['high']}  medium={tiers['medium']} ({tiers['medium']/n*100:.1f}%)  low={tiers['low']} ({tiers['low']/n*100:.1f}%)")

with open(V5 / "08_benchmark/tier_distribution_v5_noML.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species", "n", "high", "medium", "low", "medium_pct", "low_pct"])
    w.writerows(summary)
print(f"\n✅ 18 物种重聚合完成 → tier_distribution_v5_noML.tsv")
