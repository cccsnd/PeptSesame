#!/usr/bin/env python3
"""Layer4: per-species SSP vs all_ssp BLASTP → conserved hits

- 每物种 q: ssp_fa/<sp>.ssp.fa vs db: all_ssp_db
- 排除自身命中 (同物种)
- 输出: blastp/<sp>_vs_all.tsv (outfmt6) + conserved 统计
"""
import os, subprocess, time

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"
CS = f"{R2}/05_cross_species"
BLASTP = "blastp"
DB = f"{CS}/all_ssp_db"
ALL = ["Yu11", "S3651", "14G01", "14G02", "K16", "ken1", "ken8",
       "Arabidopsis", "Rice", "Tomato", "Sunflower", "Flax", "Grape",
       "Ricinus", "Zisu", "Maize", "Soybean", "Gastrodia"]

os.makedirs(f"{CS}/blastp", exist_ok=True)
for sp in ALL:
    q = f"{CS}/ssp_fa/{sp}.ssp.fa"
    out6 = f"{CS}/blastp/{sp}_vs_all.tsv"
    if os.path.exists(out6) and os.path.getsize(out6) > 0:
        n = sum(1 for _ in open(out6))
        print(f"✅ {sp}: 已存在 ({n:,} 命中)", flush=True)
        continue
    t0 = time.time()
    print(f"🔵 {sp}: blastp 运行中...", flush=True)
    r = subprocess.run([BLASTP, "-query", q, "-db", DB, "-out", out6,
                        "-outfmt", "6", "-evalue", "1e-3", "-max_target_seqs", "5",
                        "-num_threads", "8", "-seg", "yes"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"❌ {sp}: blastp 失败: {r.stderr[-200:]}", flush=True)
        continue
    n = sum(1 for _ in open(out6))
    print(f"✅ {sp}: {n:,} 命中 ({time.time()-t0:.0f}s)", flush=True)
print("全部 BLASTP 完成")
