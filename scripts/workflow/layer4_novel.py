#!/usr/bin/env python3
"""Layer4 step 3: novelty screening (vs Yu11 annotated proteins) + core set

步骤:
1. 每物种 SSP vs Yu11 注释蛋白 blastp → 有同源者排除 (novel = 无同源)
2. 核心集 = conserved (有跨物种命中) × novel (无 Yu11 同源)
输出:
- blastp/<sp>_vs_yu11.tsv (blastp outfmt6)
- conserved_novel_stats.tsv (每物种统计)
"""
import os, subprocess, time, csv

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"
CS = f"{R2}/05_cross_species"
BLASTP = "blastp"
MAKEDB = "makeblastdb"
YU11_PEP = os.environ.get("PEPTSESAME_YU11_PEP", f"{ROOT}/data/Yu11_t2t.longest_final.pep.fa")
ALL = ["Yu11", "S3651", "14G01", "14G02", "K16", "ken1", "ken8",
       "Arabidopsis", "Rice", "Tomato", "Sunflower", "Flax", "Grape",
       "Ricinus", "Zisu", "Maize", "Soybean", "Gastrodia"]

os.makedirs(f"{CS}/blastp", exist_ok=True)

# 1. Yu11 pep DB
yu11_db = f"{CS}/yu11_pep_db"
if not os.path.exists(yu11_db + ".phr"):
    print("makeblastdb Yu11 pep...")
    subprocess.run([MAKEDB, "-in", YU11_PEP, "-dbtype", "prot", "-out", yu11_db],
                   check=True, capture_output=True)
print("Yu11 pep DB 就绪")

# 2. 每物种 blastp vs Yu11 pep
stats = []
for sp in ALL:
    q = f"{CS}/ssp_fa/{sp}.ssp.fa"
    out6 = f"{CS}/blastp/{sp}_vs_yu11.tsv"
    if not (os.path.exists(out6) and os.path.getsize(out6) > 0):
        t0 = time.time()
        print(f"🔵 {sp}: blastp vs Yu11 pep...", flush=True)
        r = subprocess.run([BLASTP, "-query", q, "-db", yu11_db, "-out", out6,
                            "-outfmt", "6", "-evalue", "1e-3", "-max_target_seqs", "5",
                            "-num_threads", "8"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"❌ {sp}: {r.stderr[-200:]}", flush=True)
            continue
        print(f"✅ {sp}: 完成 ({time.time()-t0:.0f}s)", flush=True)
    # 统计
    hit_ids = set()
    with open(out6) as f:
        for line in f:
            hit_ids.add(line.split("\t")[0])
    n_ssp = sum(1 for l in open(q) if l.startswith(">"))
    n_novel = n_ssp - len(hit_ids)
    stats.append({"species": sp, "ssp": n_ssp, "yu11_hit": len(hit_ids), "novel": n_novel})
    print(f"  {sp}: SSP={n_ssp:,} Yu11同源={len(hit_ids):,} novel={n_novel:,}", flush=True)

# 3. 写统计
with open(f"{CS}/novel_stats.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=["species", "ssp", "yu11_hit", "novel"])
    w.writeheader()
    w.writerows(stats)
print(f"\n✅ novel 统计: {CS}/novel_stats.tsv")
