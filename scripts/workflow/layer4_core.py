#!/usr/bin/env python3
"""Layer4 core set: build species-prefix library from all blastp subjects

比从 ssp_fa 采样更完整: 直接收集所有 blastp 输出中出现的 subject ID,
按物种文件归属 (query 是单物种文件, 但 subject 跨物种)。
更稳的方法: 用每个物种 ssp_fa 的全部 ID (或更大采样) 构建前缀。
"""
import os, csv
from collections import Counter

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"
CS = f"{R2}/05_cross_species"
ALL = ["Yu11", "S3651", "14G01", "14G02", "K16", "ken1", "ken8",
       "Arabidopsis", "Rice", "Tomato", "Sunflower", "Flax", "Grape",
       "Ricinus", "Zisu", "Maize", "Soybean", "Gastrodia"]

# 用全部 ssp_fa ID 构建前缀 (2-token 或 # 前)
def get_prefixes_full(sp):
    fa = f"{CS}/ssp_fa/{sp}.ssp.fa"
    parts = set()
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                if "#" in sid:
                    parts.add(sid.split("#")[0])
                else:
                    toks = sid.split("_")
                    parts.add("_".join(toks[:2]))
    return parts

print("构建全量前缀库...")
SP_PREFIX = {sp: get_prefixes_full(sp) for sp in ALL}
for sp, p in SP_PREFIX.items():
    print(f"  {sp}: {len(p)} 前缀 (例: {sorted(p)[:3]})")

def species_of(sid):
    if "#" in sid:
        key = sid.split("#")[0]
    else:
        toks = sid.split("_")
        key = "_".join(toks[:2])
    for sp, prefixes in SP_PREFIX.items():
        if key in prefixes:
            return sp
    # 单 token 退路
    key1 = key.split("_")[0]
    hits = [sp for sp, prefixes in SP_PREFIX.items() if key1 in prefixes]
    return hits[0] if len(hits) == 1 else None

# 测试
tests = ["NC_050096.1_1196130_1196181_+_f0", "NC_015438.3_732870_732999_+_f0",
         "CM007890.2_75915_75957_+_f0", "1_36552_36582_+_f0",
         "NW_003018642.1_509_758_-_f2", "GWHBHOU00000001_378699_378819_+_f0"]
print("\n物种判定测试:")
for t in tests:
    print(f"  {t[:40]} → {species_of(t)}")

core_rows = []
for sp in ALL:
    cons = f"{CS}/blastp/{sp}_vs_all.tsv"
    novel = f"{CS}/blastp/{sp}_vs_yu11.tsv"
    if not (os.path.exists(cons) and os.path.exists(novel)):
        continue
    hit_yu11 = set()
    with open(novel) as f:
        for line in f:
            hit_yu11.add(line.split("\t")[0])
    best = {}
    self_hits = 0
    unknown = 0
    with open(cons) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            q, s, pident, qcovs, evalue = p[0], p[1], float(p[2]), float(p[8]), float(p[10])
            ssp = species_of(s)
            if ssp == sp:
                self_hits += 1
                continue
            if ssp is None:
                unknown += 1
                continue
            if q not in best or evalue < best[q][3]:
                best[q] = (s, pident, qcovs, evalue, ssp)
    n_core = 0
    for q, (s, pident, qcovs, evalue, ssp) in best.items():
        if q not in hit_yu11:
            # conservation_score consistent with the earlier method (pident × qcovs × e_pen × 1.5)
            e_pen = 1.0 if evalue < 1e-10 else (0.8 if evalue < 1e-5 else 0.5)
            score = round(min(1.0, (pident/100) * (qcovs/100) * e_pen * 1.5), 4)
            core_rows.append({
                "sorf_id": f"{q}|{sp}", "query_species": sp,
                "best_hit_species": ssp, "best_hit_sorf": s,
                "pident": pident, "qcovs": qcovs, "evalue": evalue,
                "conservation_score": score
            })
            n_core += 1
    print(f"{sp}: 自身={self_hits:,} 未知={unknown:,} 核心={n_core:,}")

out = f"{CS}/novel_conserved_ssp.tsv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t",
                       fieldnames=["sorf_id", "query_species", "best_hit_species",
                                   "best_hit_sorf", "pident", "qcovs", "evalue",
                                   "conservation_score"])
    w.writeheader()
    w.writerows(core_rows)
print(f"\n✅ 核心集: {out} ({len(core_rows):,} 条)")
by_sp = Counter(r["query_species"] for r in core_rows)
print("按物种:", dict(by_sp))
