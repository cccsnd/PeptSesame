#!/usr/bin/env python3
"""Layer4 cross-species conservation analysis (BLASTP)

步骤:
1. 收集 18 物种 SSP 候选 → all_ssp.fa (BLASTP 数据库)
2. 每物种 SSP vs all_ssp → conserved (cross-species hit)
3. 每物种 SSP vs Yu11 注释蛋白 → novel (无同源)
4. 核心集 = conserved × novel

依赖: blastp (需在 PATH 或指定)
输出: results_v2/05_cross_species/
"""
import os, sys, subprocess, json

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"
CS = f"{R2}/05_cross_species"
os.makedirs(CS, exist_ok=True)

with open(f"{R2}/species_manifest.json") as f:
    MANIFEST = json.load(f)

def get_ssp_fa(sp):
    """提取 SSP 候选序列 → fasta (供 blastp)"""
    cls = f"{R2}/03_layer3_classify/{sp}/classified_sorfs.tsv"
    fa = f"{R2}/01_layer1_sixframe/{sp}/sorfs.fa"
    out = f"{CS}/ssp_fa/{sp}.ssp.fa"
    os.makedirs(f"{CS}/ssp_fa", exist_ok=True)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    # 读 SSP id
    ssp_ids = set()
    with open(cls) as f:
        import csv
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] == "True":
                ssp_ids.add(r["seq_id"])
    # 从 FASTA 提取
    seqs = {}
    cur = None; buf = []
    with open(fa) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur: seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if cur: seqs[cur] = "".join(buf)
    n = 0
    with open(out, "w") as f:
        for sid in seqs:
            key = sid if sid in ssp_ids else None
            if key:
                f.write(f">{key}\n{seqs[key]}\n")
                n += 1
    print(f"  {sp}: SSP FASTA {n} 条 → {out}")
    return out

if __name__ == "__main__":
    # 只收集 SSP (不做 blastp, 等 Layer2/3 全部完成)
    for sp in MANIFEST:
        cls = f"{R2}/03_layer3_classify/{sp}/classified_sorfs.tsv"
        if os.path.exists(cls):
            print(f"=== {sp} ===")
            get_ssp_fa(sp)
        else:
            print(f"⚠️ {sp} Layer3 未完成, 跳过")
