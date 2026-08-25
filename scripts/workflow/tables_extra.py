#!/usr/bin/env python3
"""tables_extra.py — 补充表 S2 (SSP 序列) + S4 (Top100) 生成 (current pipeline)

流式实现: 先收集 is_ssp=True 的 id 集合 (小), 再流式扫 FASTA 只输出命中,
避免全量序列载入内存。S4 = Yu11 评分 Top100 (aggregated_score 降序 + 家族)。
"""
import csv
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results"
OUT = R2 / "09_tables"
OUT.mkdir(parents=True, exist_ok=True)

FAM_COUNTS = OUT / "family_counts.tsv"
SPECIES = [r["species"] for r in csv.DictReader(open(FAM_COUNTS), delimiter="\t")]


def main():
    n_total = 0
    with open(OUT / "TableS2_SSP_sequences.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Genome", "sORF_ID", "SSP_Family", "Length_aa",
                    "Peptide_Sequence"])
        for sp in SPECIES:
            cls = R2 / f"03_layer3_classify/{sp}/classified_sorfs.tsv"
            fa = R2 / f"01_layer1_sixframe/{sp}/sorfs.fa"
            if not (cls.exists() and fa.exists()):
                print(f"⚠️ {sp} 输入缺失")
                continue
            # 1) 收集 is_ssp=True 的 id + 家族 (小集合)
            keep = {}
            with open(cls) as fsrc:
                for r in csv.DictReader(fsrc, delimiter="\t"):
                    if r["is_ssp"] == "True":
                        keep[r["seq_id"]] = r["ssp_families"] or "-"
                        keep[r["seq_id"].replace("sORF_", "")] = keep[r["seq_id"]]
            # 2) 流式扫 FASTA, 只输出命中
            n = 0
            sid = None
            buf = []
            with open(fa) as ffa:
                for line in ffa:
                    line = line.strip()
                    if line.startswith(">"):
                        if sid and sid in keep:
                            seq = "".join(buf)
                            w.writerow([sp, sid, keep[sid], len(seq), seq])
                            n += 1
                        sid = line[1:].split()[0]
                        buf = []
                    else:
                        buf.append(line)
                if sid and sid in keep:
                    seq = "".join(buf)
                    w.writerow([sp, sid, keep[sid], len(seq), seq])
                    n += 1
            print(f"  {sp}: {n:,} SSP 序列", flush=True)
            n_total += n
    print(f"✅ TableS2_SSP_sequences.tsv (总计 {n_total:,} 条)")

    # S4: Yu11 Top100
    scored = R2 / "02_layer2_scoring/Yu11/scored_sorfs.tsv"
    cls = R2 / "03_layer3_classify/Yu11/classified_sorfs.tsv"
    fam_of = {}
    with open(cls) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fam_of[r["seq_id"]] = r["ssp_families"] or "-"
    top = []
    with open(scored) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            top.append((float(r["aggregated_score"]), r["seq_id"],
                        r["confidence"]))
    top.sort(reverse=True)
    with open(OUT / "TableS4_top100.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Rank", "sORF_ID", "SSP_Family", "Aggregated_Score",
                    "Confidence"])
        for i, (score, sid, conf) in enumerate(top[:100], 1):
            w.writerow([i, sid, fam_of.get(sid, "-"), round(score, 4), conf])
    print(f"✅ TableS4_top100.tsv (100 行)")


if __name__ == "__main__":
    main()
