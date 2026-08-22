import os
#!/usr/bin/env python3
"""enrichment_matched.py — 匹配背景富集 (Supplementary Table S9)

Background for core-set enrichment must be matched by length/composition/family.
方法: 分层 permutation — 按 (家族 × 长度 bin 10aa) 分层, 从非 DE 的
Yu11 SSP 中按 DE 集合的层比例抽取 240 个背景 × 10,000 次, 统计核心集
命中数分布, 观测值 55 的置换 p 值 + 期望命中率。
输出: results_v2/08_benchmark/enrichment_matched.tsv
"""
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results_v2"
OUT = R2 / "08_benchmark"
N_PERM = 10000
RNG = random.Random(42)


def norm_id(s):
    return s.replace("sORF_", "").split("|")[0]


def main():
    # 1. 核心集
    core = set()
    with open(R2 / "05_cross_species/novel_conserved_ssp.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            core.add(norm_id(r["sorf_id"]))

    # 2. DE sORF
    de_genes = set()
    with open(R2 / "06_expression/de_candidates_strict.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            de_genes.add(r["gene"])
    de_sorfs = set()
    with open(R2 / "06_expression/sorf_gene_map.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["gene"] in de_genes:
                de_sorfs.add(norm_id(r["sorf_id"]))
    print(f"DE sORF: {len(de_sorfs)}")

    # 3. Yu11 全 SSP: 家族 (首位) + 长度
    rows = []
    with open(R2 / "03_layer3_classify/Yu11/classified_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] != "True":
                continue
            sid = norm_id(r["seq_id"])
            fam = (r["ssp_families"] or "?").split(";")[0]
            ln = int(r["length_aa"])
            lbin = (ln // 10) * 10
            rows.append({
                "id": sid, "fam": fam, "lbin": lbin,
                "de": sid in de_sorfs, "core": sid in core,
            })
    n_all = len(rows)
    n_de = sum(1 for r in rows if r["de"])
    n_de_core = sum(1 for r in rows if r["de"] and r["core"])
    print(f"Yu11 SSP: {n_all}, DE: {n_de}, DE∩核心集: {n_de_core}")

    # 4. 分层: (fam, lbin) -> [非 DE 成员]; DE 集合的层分布
    strata = defaultdict(list)
    de_strata = Counter()
    for r in rows:
        key = (r["fam"], r["lbin"])
        if r["de"]:
            de_strata[key] += 1
        else:
            strata[key].append(r["id"])
    # 检查每层非 DE 池够不够 (不足的层合并到最近层: 简单做法: 不足则从全局池补)
    total_pool = [r["id"] for r in rows if not r["de"]]
    # 5. 置换
    obs = n_de_core
    hits = 0
    n_under = 0
    exp_sum = 0.0
    for it in range(N_PERM):
        picked = []
        for key, k in de_strata.items():
            pool = strata.get(key, [])
            if len(pool) < k:
                pool = pool + total_pool  # 补充 (极少触发)
            picked.extend(RNG.sample(pool, k))
        h = sum(1 for sid in picked if sid in core)
        exp_sum += h / len(picked)
        if h >= obs:
            hits += 1
    p = (hits + 1) / (N_PERM + 1)
    exp_rate = exp_sum / N_PERM
    obs_rate = obs / n_de
    print(f"观测: {obs}/{n_de} = {obs_rate*100:.1f}%")
    print(f"置换期望: {exp_rate*100:.2f}% (n={N_PERM})")
    print(f"置换 p (≥{obs}): {p:.4f}")

    with open(OUT / "enrichment_matched.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_de_sorfs\t{n_de}\n")
        f.write(f"observed_core_hits\t{obs}\n")
        f.write(f"observed_rate\t{obs_rate:.4f}\n")
        f.write(f"permutation_n\t{N_PERM}\n")
        f.write(f"expected_rate\t{exp_rate:.4f}\n")
        f.write(f"permutation_p\t{p:.4f}\n")
        f.write(f"note\tstratified by family x 10aa-length-bin; "
                f"permutation p = P(hits >= obs)\n")
    print(f"✅ {OUT / 'enrichment_matched.tsv'}")


if __name__ == "__main__":
    main()
