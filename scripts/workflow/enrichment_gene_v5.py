import os
#!/usr/bin/env python3
"""D. 基因级分层置换富集 — R5/R6 要求 (以 host gene 为统计单位, 控制非独立性)

背景: 240 records = 240 unique sORF (零多关联), 但 11 个基因含多候选 →
      association-record 级置换低估基因聚类非独立性。
方法: 以 221 DE 基因为观测集, 从非 DE 的 Yu11 motif-compatible 宿主基因中
      按 (家族 × 10-aa 长度 bin) 分层抽样 221 个背景基因 × 10,000 次,
      统计核心集命中基因数, 观测 54/221 的置换 P 值。
输出: results/08_benchmark/enrichment_gene_v20260825.tsv
"""
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

V5 = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "results"
OUT = V5 / "08_benchmark/enrichment_gene_v20260825.tsv"
N_PERM = 10000
RNG = random.Random(42)


def norm_id(s):
    return s.replace("sORF_", "").split("|")[0]


def main():
    # 核心集
    core = set()
    with open(V5 / "05_cross_species/novel_conserved_ssp.tsv") as f:
        next(f)
        for line in f:
            core.add(line.split("\t")[0].split("|")[0].replace("sORF_", ""))
    print(f"核心集: {len(core):,}")

    # DE 基因
    de_genes = set()
    with open(V5 / "06_expression/de_candidates_strict.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            de_genes.add(r["gene"])
    print(f"DE 基因: {len(de_genes)}")

    # 基因 → (family, length) 层 (取该基因候选的首位家族 + 最长 sORF 长度)
    gene_layer = {}          # gene -> (fam, lbin)
    gene_core = {}           # gene -> 是否核心集命中
    with open(V5 / "03_layer3_classify/Yu11/classified_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] != "True":
                continue
            sid = norm_id(r["seq_id"])
            fam = (r["ssp_families"] or "?").split(";")[0]
            ln = int(r["length_aa"])
            lbin = (ln // 10) * 10
            gene = None
            # 从 sorf_gene_map 找宿主基因 (仅 Yu11 候选)
            # (预加载 map 加速)
            gene_layer.setdefault(sid, (fam, lbin))
    # sorf_gene_map: sorf -> gene
    sorf2gene = {}
    with open(V5 / "06_expression/sorf_gene_map.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sorf2gene.setdefault(norm_id(r["sorf_id"]), r["gene"])

    de_layer = {}
    de_core_hit = 0
    de_genes_with_core = set()
    for sid, (fam, lbin) in gene_layer.items():
        g = sorf2gene.get(sid)
        if g is None:
            continue
        if g in de_genes:
            de_layer[g] = (fam, lbin)
            if sid in core:
                de_core_hit += 1
                de_genes_with_core.add(g)
    print(f"DE 基因层: {len(de_layer)}, 核心集命中基因: {len(de_genes_with_core)} ({len(de_genes_with_core)/len(de_layer)*100:.1f}%)")

    # 背景池: 非 DE 基因 (全部基因的层)
    all_genes = set(sorf2gene.values())
    non_de = all_genes - de_genes
    pool = defaultdict(list)
    for g in non_de:
        # 该基因的层 (取第一个 sORF 的层; 多 sORF 基因各层都加入池以保持层比例)
        for sid, (fam, lbin) in gene_layer.items():
            if sorf2gene.get(sid) == g:
                pool[(fam, lbin)].append(g)
                break
    # 观测比例按 DE 层分布抽样
    de_dist = Counter(de_layer.values())
    obs = len(de_genes_with_core)
    n_ge = 0
    for _ in range(N_PERM):
        bg = set()
        for (fam, lbin), cnt in de_dist.items():
            cands = pool.get((fam, lbin), [])
            if len(cands) >= cnt:
                bg.update(RNG.sample(cands, cnt))
            else:
                bg.update(cands)
        hits = sum(1 for g in bg if any(sid in core for sid, gg in sorf2gene.items() if gg == g and (sorf2gene.get(sid) == g)))
        # 简化: 预计算基因→核心
        if hits >= obs:
            n_ge += 1
    p = (1 + n_ge) / (1 + N_PERM)
    exp_rate = None
    # 期望率: 抽样背景的核心命中均值 (用独立 500 次估计)
    exp_hits = 0
    for _ in range(500):
        bg = set()
        for (fam, lbin), cnt in de_dist.items():
            cands = pool.get((fam, lbin), [])
            bg.update(RNG.sample(cands, cnt) if len(cands) >= cnt else cands)
        exp_hits += sum(1 for g in bg if gene_core_flag(g))
    exp_rate = exp_hits / 500 / len(de_genes)

    print(f"观测: {obs}/{len(de_genes)} = {obs/len(de_genes)*100:.1f}%")
    print(f"期望率 (500 次估计): {exp_rate*100:.2f}%")
    print(f"置换 P: {p:.4f}")
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerows([
            ["unit", "host gene (221)"],
            ["n_genes", len(de_genes)],
            ["observed_core_genes", obs],
            ["observed_rate", f"{obs/len(de_genes):.4f}"],
            ["expected_rate", f"{exp_rate:.4f}"],
            ["permutation_n", N_PERM],
            ["permutation_p", f"{p:.4f}"],
            ["note", "stratified by family x 10-aa length bin at host-gene level"],
        ])
    print(f"✅ {OUT}")


def gene_core_flag(g):
    return g in _gene_core_cache


_gene_core_cache = set()
if __name__ == "__main__":
    # 预计算基因→核心 (避免循环内重复扫描)
    V5 = Path("{V5}")
    core = set()
    with open(V5 / "05_cross_species/novel_conserved_ssp.tsv") as f:
        next(f)
        for line in f:
            core.add(line.split("\t")[0].split("|")[0].replace("sORF_", ""))
    sorf2gene = {}
    with open(V5 / "06_expression/sorf_gene_map.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sorf2gene.setdefault(r["sorf_id"].replace("sORF_", "").split("|")[0], r["gene"])
    for sid, g in sorf2gene.items():
        if sid in core:
            _gene_core_cache.add(g)
    main()
