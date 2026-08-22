#!/usr/bin/env python3
"""Layer4: extract per-species SSP FASTA (background, per-species to avoid timeout)"""
import os, csv

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"
CS = f"{R2}/05_cross_species"
ALL = ["Yu11", "S3651", "14G01", "14G02", "K16", "ken1", "ken8",
       "Arabidopsis", "Rice", "Tomato", "Sunflower", "Flax", "Grape",
       "Ricinus", "Zisu", "Maize", "Soybean", "Gastrodia"]
os.makedirs(f"{CS}/ssp_fa", exist_ok=True)

for sp in ALL:
    cls = f"{R2}/03_layer3_classify/{sp}/classified_sorfs.tsv"
    fa = f"{R2}/01_layer1_sixframe/{sp}/sorfs.fa"
    out = f"{CS}/ssp_fa/{sp}.ssp.fa"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print(f"✅ {sp}: 已存在", flush=True)
        continue
    # 读 SSP id
    ssp_ids = set()
    with open(cls) as f:
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
        for sid, seq in seqs.items():
            # Layer3 seq_id 带 sORF_ 前缀, FASTA id 无前缀 — 两种都匹配
            if sid in ssp_ids or f"sORF_{sid}" in ssp_ids:
                f.write(f">{sid}\n{seq}\n")
                n += 1
    print(f"✅ {sp}: {n:,} SSP → {out}", flush=True)
print("全部提取完成")
