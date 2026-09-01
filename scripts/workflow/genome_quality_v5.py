#!/usr/bin/env python3
"""G. 18 基因组装配质量表 (可计算指标) — R5/R6 要求

从各基因组 FASTA 计算: 总长 / 序列数 / N50 / 最长序列
(软掩蔽状态从文件名与来源记录标注; BUSCO 未运行 — 明确标注)
输出: results/09_tables/TableS14_genome_quality.tsv
"""
import os
import csv
from pathlib import Path

GENOMES = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "results/00_inputs/genomes"
OUT = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "results/09_tables/TableS14_genome_quality.tsv"


def stats(fa):
    n = 0
    total = 0
    lens = []
    cur = 0
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                    total += cur
                    n += 1
                    cur = 0
                n += 0  # 序列计数在 close
            else:
                cur += len(line.strip())
        if cur:
            lens.append(cur)
            total += cur
            n += 1
    lens.sort(reverse=True)
    half = total / 2
    acc = 0
    n50 = 0
    for L in lens:
        acc += L
        if acc >= half:
            n50 = L
            break
    return total, n, n50, lens[0] if lens else 0


rows = []
for fa in sorted(GENOMES.glob("*.fasta")):
    sp = fa.stem
    if sp in ("Maize_B73", "Yuzhi11"):  # 别名文件 (Maize/Yu11 副本), 不计入 18 物种
        continue
    total, n, n50, longest = stats(fa)
    rows.append([sp, f"{total/1e6:.2f}", n, f"{n50/1e6:.2f}", f"{longest/1e6:.2f}", "unmasked (scanned as-is)"])
    print(f"{sp:<14} {total/1e6:8.2f} Mb  {n:5d} 条  N50={n50/1e6:7.2f} Mb  最长={longest/1e6:6.2f} Mb")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species", "size_Mb", "n_sequences", "N50_Mb", "longest_Mb", "masking_note"])
    w.writerows(rows)
print(f"\n✅ {OUT}")
print("注: BUSCO/注释完整性未运行 (需各基因组 lineage DB); 掩蔽状态 = 扫描输入即原样 (软掩蔽小写转大写, 见 Methods)")
