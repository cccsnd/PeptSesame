#!/usr/bin/env python3
"""v2 全流程重跑 — Layer1 六框翻译 (从原始基因组 FASTA, 输出到 results_v2)

严格从零: 输入 = 原始基因组 FASTA + GFF (只读); 输出 = results_v2/01_layer1_sixframe/
不复用 results/ 任何文件

用法: python run_layer1.py <species> <genome_fasta> [gff]
"""
import os, sys, time, logging

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
sys.path.insert(0, f"{ROOT}/pipeline/layer1_sixframe")
from sixframe import SixFrameTranslator

def run(species, genome_fasta, gff=None, min_len=30, max_len=300):
    t0 = time.time()
    out_dir = f"{ROOT}/results_v2/01_layer1_sixframe/{species}"
    os.makedirs(out_dir, exist_ok=True)
    bed = f"{out_dir}/sorfs.bed"
    fa = f"{out_dir}/sorfs.fa"

    translator = SixFrameTranslator(
        genome_fasta=genome_fasta,
        cds_gff=gff,
        min_orf_len=min_len,
        max_orf_len=max_len,
        out_bed=bed,
        out_pep_fasta=fa,
    )
    translator.run()
    n = sum(1 for _ in open(bed)) - 1 if os.path.exists(bed) else 0
    print(f"[{species}] Layer1 完成: {n:,} sORF ({time.time()-t0:.0f}s)", flush=True)
    return n

if __name__ == "__main__":
    species = sys.argv[1]
    genome = sys.argv[2]
    gff = sys.argv[3] if len(sys.argv) > 3 else None
    run(species, genome, gff)
