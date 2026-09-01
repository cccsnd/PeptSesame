import os
#!/usr/bin/env python3
"""重算 classified_sorfs.tsv 的 score/tier (去 ML, 5 通道归一化) + SSP 候选 tier 分布"""
import csv
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
V5 = ROOT / "results"
SPECIES = ["Yu11", "14G01", "14G02", "K16", "ken1", "ken8", "S3651",
           "Arabidopsis", "Rice", "Maize", "Soybean", "Tomato", "Grape",
           "Ricinus", "Flax", "Sunflower", "Zisu", "Gastrodia"]

W = {"sequence_features": 0.20/0.75, "conservation": 0.15/0.75,
     "expression": 0.10/0.75, "structural": 0.15/0.75, "motif": 0.15/0.75}


def tier_of(s):
    return "high" if s >= 0.7 else ("medium" if s >= 0.4 else "low")


out_rows = []
for sp in SPECIES:
    p = V5 / f"03_layer3_classify/{sp}/classified_sorfs.tsv"
    if not p.exists():
        continue
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    ssp_tiers = {"high": 0, "medium": 0, "low": 0}
    ssp_n = 0
    for r in rows:
        s = sum(W[k] * float(r[k]) for k in W)
        r["aggregated_score"] = f"{s:.6f}"
        r["confidence"] = tier_of(s)
        if r["is_ssp"] == "True":
            ssp_n += 1
            ssp_tiers[tier_of(s)] += 1
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    if ssp_n:
        out_rows.append((sp, ssp_n, ssp_tiers["high"], ssp_tiers["medium"], ssp_tiers["low"],
                         f"{ssp_tiers['medium']/ssp_n*100:.1f}%", f"{ssp_tiers['low']/ssp_n*100:.1f}%"))
        print(f"  {sp}: SSP n={ssp_n}  high={ssp_tiers['high']}  medium={ssp_tiers['medium']} ({ssp_tiers['medium']/ssp_n*100:.1f}%)  low={ssp_tiers['low']} ({ssp_tiers['low']/ssp_n*100:.1f}%)")

with open(V5 / "08_benchmark/ssp_tier_distribution_v5_noML.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species", "n_ssp", "high", "medium", "low", "medium_pct", "low_pct"])
    w.writerows(out_rows)
print(f"\n✅ SSP tier 分布 → ssp_tier_distribution_v5_noML.tsv")
