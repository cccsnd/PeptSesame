#!/usr/bin/env python3
"""Layer2 scoring (per-species strict pipeline)

用法: python run_layer2.py <species> [--max-sorfs N]
输出: results/02_layer2_scoring/<species>/scored_sorfs.tsv
自检: 行数 = Layer1 BED 行数 (除 header); 分数范围 0-1; 无空值
"""
import os, sys, csv, time, argparse

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results"

sys.path.insert(0, ROOT)  # 以包方式导入 pipeline.layer2_scoring
from pipeline.layer2_scoring.scoring_core import EvidenceScorer

def main():
    p = argparse.ArgumentParser()
    p.add_argument("species")
    p.add_argument("--max-sorfs", type=int, default=0)
    args = p.parse_args()
    sp = args.species

    d = f"{R2}/01_layer1_sixframe/{sp}"
    bed = f"{d}/sorfs.bed"
    fa = f"{d}/sorfs.fa"
    if not (os.path.exists(bed) and os.path.exists(fa)):
        print(f"❌ {sp} Layer1 产物缺失: {bed} / {fa}")
        sys.exit(1)

    # 预期行数 (BED 总行数 - header)
    n_bed = sum(1 for _ in open(bed)) - 1
    print(f"=== {sp} Layer2 评分 ===")
    print(f"Layer1 BED 行数: {n_bed:,}")

    from Bio import SeqIO
    pep_map = {}
    for r in SeqIO.parse(fa, "fasta"):
        pep_map[r.id] = str(r.seq)
    print(f"FASTA 序列: {len(pep_map):,}")

    # 读 BED (匹配 FASTA)
    sorfs = []
    with open(bed) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            name = parts[3]
            seq = pep_map.get(name, "") or pep_map.get(name.replace("sORF_", ""), "")
            if not seq:
                continue
            sorfs.append({
                "seq_id": name, "sequence": seq,
                "chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                "strand": parts[5], "length_nt": int(parts[7]), "length_aa": int(parts[8]),
            })
            if args.max_sorfs and len(sorfs) >= args.max_sorfs:
                break
    print(f"加载 sORF: {len(sorfs):,} (匹配 FASTA)")
    if args.max_sorfs == 0 and len(sorfs) != n_bed:
        print(f"⚠️ 加载数 {len(sorfs):,} != BED 行数 {n_bed:,} — 有未匹配!")

    # 评分 (rule-only, 无 ML — 与稿件口径一致)
    t0 = time.time()
    scorer = EvidenceScorer(ml_model_path=None)
    out_dir = f"{R2}/02_layer2_scoring/{sp}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/scored_sorfs.tsv"
    n = len(sorfs)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["seq_id", "chrom", "start", "end", "strand", "length_nt", "length_aa",
                    "sequence_features", "conservation", "expression", "structural",
                    "motif", "ml_score", "aggregated_score", "confidence"])
        for i, s in enumerate(sorfs):
            r = scorer.evaluate(s["seq_id"], s["sequence"])
            w.writerow([s["seq_id"], s["chrom"], s["start"], s["end"], s["strand"],
                        s["length_nt"], s["length_aa"],
                        r["sub_scores"]["sequence_features"]["combined"],
                        r["sub_scores"]["conservation"]["combined"],
                        r["sub_scores"]["expression"]["combined"],
                        r["sub_scores"]["structural"]["combined"],
                        r["sub_scores"]["motif"]["combined"],
                        r["sub_scores"]["ml"]["combined"],
                        r["aggregated_score"], r["confidence"]])
            if (i + 1) % 50000 == 0:
                print(f"  {i+1:,}/{n:,} ({(i+1)/(time.time()-t0):.0f}/s)", flush=True)
    print(f"评分完成: {n:,} → {out_path} ({time.time()-t0:.0f}s)")

    # 自检
    n_out = sum(1 for _ in open(out_path)) - 1
    assert n_out == n, f"输出行数 {n_out} != {n}"
    print(f"✅ {sp} Layer2 完成: {n_out:,} 行, 自检通过")

if __name__ == "__main__":
    main()
