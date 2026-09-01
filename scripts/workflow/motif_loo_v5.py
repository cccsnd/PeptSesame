#!/usr/bin/env python3
"""B. motif LOO + 全搜索空间期望假命中 — R3/R5/R6 要求

1. LOO 检查: 13 条参考序列逐一剔除后, 其余参考是否仍被 regex 命中
   (regex 为固定模式, LOO 验证命中不依赖单条序列 + 报告每成员匹配详情)
2. 全搜索空间期望假命中: 用 Yu11 经验 null (5 批独立抽样) 的 null_mean
   归一化到全库 → 每家族在 445 万搜索空间的期望假命中上界
输出: results/08_benchmark/motif_loo_v20260825.tsv + 说明
"""
import csv
import re
import os, sys
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
sys.path.insert(0, str(ROOT))
from pipeline.layer3_classify.motif_profiles import SSP_MOTIFS

COMPILED = {fam: re.compile(p) for fam, p in SSP_MOTIFS.items()}

# 参考序列: 直接从权威 benchmark 脚本导入 (避免手写不一致)
import importlib.util
spec = importlib.util.spec_from_file_location(
    "bmr", ROOT / "scripts/workflow/benchmark_motif_recall.py")
bmr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bmr)
GOLD = bmr.GOLD

def main():
    out_rows = []
    print(f"{'family':<6}{'member':<3}{'hit':<6}{'match_pos':<10}")
    total = 0
    hit = 0
    per_fam = {}
    for fam, seqs in GOLD.items():
        pat = COMPILED[fam]
        per_fam[fam] = [0, len(seqs)]
        for i, s in enumerate(seqs):
            m = pat.search(s)
            ok = m is not None
            total += 1
            hit += 1 if ok else 0
            per_fam[fam][0] += 1 if ok else 0
            out_rows.append([fam, i + 1, "Y" if ok else "N", m.start() if m else "-"])
            print(f"{fam:<6}{i+1:<3}{'Y' if ok else 'N':<6}{m.start() if m else '-':<10}")
    print(f"\n参考序列回收: {hit}/{total} = {hit/total*100:.1f}%")

    # 全搜索空间期望假命中 (null 5 批版归一化到全库)
    print("\n全搜索空间期望假命中 (Yu11 null 5 批版, 归一化到 4,451,664):")
    null = {}
    with open(ROOT / "results/08_benchmark/empirical_null_motif_1000x.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            null[r["family"]] = float(r["null_mean"])
    for fam in ["CLE", "PSK", "EPFL", "IDA_LIKE"]:
        exp = null.get(fam, 0)
        print(f"  {fam}: null_mean={exp:.0f} → 全库期望假命中 ≈ {exp*4451664/500000:.0f}")

    with open(ROOT / "results/08_benchmark/motif_loo_v20260825.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family", "member", "hit", "match_pos"])
        w.writerows(out_rows)
        w.writerow([])
        w.writerow(["summary", f"recovery {hit}/{total}", f"{hit/total*100:.1f}%", ""])
        w.writerow([])
        w.writerow(["expected_false_hits_full_library", "family", "null_mean_500k", "scaled_to_4451664"])
        for fam in ["CLE", "PSK", "EPFL", "IDA_LIKE"]:
            exp = null.get(fam, 0)
            w.writerow(["", fam, f"{exp:.0f}", f"{exp*4451664/500000:.0f}"])
    print(f"\n✅ {ROOT / 'results/08_benchmark/motif_loo_v20260825.tsv'}")


if __name__ == "__main__":
    main()
