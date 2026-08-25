#!/usr/bin/env python3
"""Layer4 step 2: BLASTP conservation + novelty screening → core set

用法: python run_layer4.py [--skip-blastp]
依赖: blastp + makeblastdb (graphpan env)

输出:
- results/05_cross_species/conserved_hits.tsv (每 SSP 的跨物种命中)
- results/05_cross_species/novel_vs_yu11.tsv (vs Yu11 注释蛋白)
- results/05_cross_species/novel_conserved_ssp.tsv (核心集)
"""
import os, sys, subprocess, csv, json

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results"
CS = f"{R2}/05_cross_species"
os.makedirs(CS, exist_ok=True)

BLASTP = "blastp"
MAKEDB = "makeblastdb"
YU11_PEP = os.environ.get("PEPTSESAME_YU11_PEP", f"{ROOT}/data/Yu11_t2t.longest_final.pep.fa")

SPECIES = ["Yu11", "S3651", "14G01", "14G02", "K16", "ken1", "ken8",
           "Arabidopsis", "Rice", "Tomato", "Sunflower", "Flax", "Grape",
           "Ricinus", "Zisu", "Maize", "Soybean", "Gastrodia"]

def extract_ssp(sp):
    cls = f"{R2}/03_layer3_classify/{sp}/classified_sorfs.tsv"
    fa = f"{R2}/01_layer1_sixframe/{sp}/sorfs.fa"
    out = f"{CS}/ssp_fa/{sp}.ssp.fa"
    os.makedirs(f"{CS}/ssp_fa", exist_ok=True)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    if not os.path.exists(cls):
        return None
    ssp_ids = set()
    with open(cls) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] == "True":
                ssp_ids.add(r["seq_id"])
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
    with open(out, "w") as f:
        for sid, seq in seqs.items():
            if sid in ssp_ids:
                f.write(f">{sid}\n{seq}\n")
    return out

def main():
    # 1. 提取 + 合并
    print("=== 步骤1: 提取 SSP ===")
    fa_list = []
    for sp in SPECIES:
        p = extract_ssp(sp)
        if p:
            fa_list.append(p)
    all_ssp = f"{CS}/all_ssp.fa"
    with open(all_ssp, "w") as fout:
        for p in fa_list:
            with open(p) as fin:
                fout.write(fin.read())
    n_all = sum(1 for l in open(all_ssp) if l.startswith(">"))
    print(f"all_ssp.fa: {n_all:,} 条")

    # 2. makeblastdb
    db_prefix = f"{CS}/all_ssp_db"
    if not os.path.exists(db_prefix + ".phr"):
        print("=== 步骤2: makeblastdb ===")
        subprocess.run([MAKEDB, "-in", all_ssp, "-dbtype", "prot", "-out", db_prefix],
                       check=True, capture_output=True)
    print("DB 就绪")

    # 3. 每物种 blastp vs all_ssp (找 conserved)
    print("=== 步骤3: BLASTP 保守性 ===")
    all_hits = []
    for sp in SPECIES:
        q = f"{CS}/ssp_fa/{sp}.ssp.fa"
        out6 = f"{CS}/blastp/{sp}_vs_all.tsv"
        if not os.path.exists(out6) or os.path.getsize(out6) == 0:
            cmd = [BLASTP, "-query", q, "-db", db_prefix, "-out", out6,
                   "-outfmt", "6", "-evalue", "1e-3", "-max_target_seqs", "5",
                   "-num_threads", "8", "-seg", "yes"]
            print(f"  {sp}: blastp...")
            subprocess.run(cmd, check=True, capture_output=True)
        n_hits = sum(1 for _ in open(out6))
        print(f"  {sp}: {n_hits:,} 命中")
        all_hits.append((sp, out6, n_hits))

    # 4. vs Yu11 注释蛋白 (novel)
    print("=== 步骤4: novel 筛选 (vs Yu11 注释蛋白) ===")
    yu11_db = f"{CS}/yu11_pep_db"
    if not os.path.exists(yu11_db + ".phr"):
        subprocess.run([MAKEDB, "-in", YU11_PEP, "-dbtype", "prot", "-out", yu11_db],
                       check=True, capture_output=True)
    novel_count = {}
    for sp, q_path, _ in all_hits:
        out6 = f"{CS}/blastp/{sp}_vs_yu11.tsv"
        if not os.path.exists(out6) or os.path.getsize(out6) == 0:
            cmd = [BLASTP, "-query", q_path, "-db", yu11_db, "-out", out6,
                   "-outfmt", "6", "-evalue", "1e-3", "-max_target_seqs", "5",
                   "-num_threads", "8"]
            print(f"  {sp}: blastp vs Yu11 pep...")
            subprocess.run(cmd, check=True, capture_output=True)
        hit_ids = set()
        with open(out6) as f:
            for line in f:
                hit_ids.add(line.split("\t")[0])
        n_ssp = sum(1 for l in open(q_path) if l.startswith(">"))
        novel_count[sp] = (n_ssp, len(hit_ids))
        print(f"  {sp}: SSP={n_ssp:,} 有同源={len(hit_ids):,} novel={n_ssp-len(hit_ids):,}")

    print("\n=== 步骤5: 核心集 (conserved × novel) ===")
    # next: parse conserved hits + novel -> intersection

if __name__ == "__main__":
    main()
