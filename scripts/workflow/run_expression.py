#!/usr/bin/env python3
"""run_expression.py — current pipeline 表达/DE 分析 (results_v2/06_expression)

Method identical to the established screen_resistance_ssp.py (current pipeline):
  1. SSP 候选 = Layer3 is_ssp=True (motif library, Yu11 = 5,289)
  2. sORF→基因: 坐标重叠, 取最大 overlap (GFF gene 特征, 染色体 yu11#1#ChrNN 直接匹配)
  3. FPKM 匹配: gid / gid.1 双尝试
  4. DE: Yu11 茎(Sm) 48h 4C vs CK (3v3 均值), log2FC = log2((t+0.1)/(c+0.1))
  5. 严格阈值: |log2FC|>1 且 max(t48,c48)>0.5; 宽松: |log2FC|>0.5 且 max(t48,c48)>0.5
  6. 家族归属: 聚合 Layer3 家族列 (is_cle..is_rgf, is_ida_like 单独标注)
  7. 核心集交集: sorf_id ∈ 05_cross_species/novel_conserved_ssp.tsv
  8. FDR: gene-level padj lookup (reusing the revised deseq2_approx_yu11_sm48.tsv,
     count 矩阵为基因级, 与 sORF 版本无关; padj<0.05 为 FDR 通过)

输出 (results_v2/06_expression/):
  - sorf_gene_map.tsv            全部 SSP sORF→基因映射 + 家族 + 核心集标记
  - expressed_sorfs.tsv          有 FPKM 数据的 SSP sORF (基因级)
  - de_candidates_strict.tsv  |log2FC|>1 候选 (主文件)
  - de_candidates_loose.tsv   |log2FC|>0.5 宽松版
  - de_classification.tsv     基因级 DE 分类 (家族归属/FDR/核心集)
  - evidence_chain.md         证据链数字 (供稿件引用)
"""
from __future__ import annotations

import csv
import os
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results_v2"
OUT = R2 / "06_expression"
OUT.mkdir(parents=True, exist_ok=True)

CLASSIFIED = R2 / "03_layer3_classify/Yu11/classified_sorfs.tsv"
GFF = R2 / "00_inputs/annotations/Yu11.gff3"
FPKM_CSV = R2 / "00_inputs/rnaseq/FPKM_matrix.csv"
CORE_SET = R2 / "05_cross_species/novel_conserved_ssp.tsv"
FDR_TABLE = ROOT / "results/paper1/fdr_analysis/deseq2_approx_yu11_sm48.tsv"

FAM_COLS = ["is_cle", "is_ralf", "is_cep", "is_psk", "is_psy1",
            "is_ida", "is_epfl", "is_rgf"]


def load_core_set() -> set:
    ids = set()
    with open(CORE_SET) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ids.add(r["sorf_id"])
    return ids


def load_ssp_rows() -> list:
    """Layer3 Yu11, is_ssp=True"""
    rows = []
    with open(CLASSIFIED) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("is_ssp") == "True":
                rows.append(r)
    return rows


def load_genes(gff_path) -> dict:
    """chrom -> sorted [(start, end, gid)]"""
    gene_by_chrom = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = parts[8]
            gid = ""
            for tok in attr.split(";"):
                if tok.startswith("ID="):
                    gid = tok[3:]
                    break
            if gid:
                gene_by_chrom[parts[0]].append(
                    (int(parts[3]), int(parts[4]), gid))
    for chrom in gene_by_chrom:
        gene_by_chrom[chrom].sort()
    return gene_by_chrom


def map_sorf_to_gene(rows, gene_by_chrom) -> dict:
    """sorf_id -> (gid, overlap_bp), 最大 overlap"""
    mapping = {}
    for r in rows:
        chrom = r["chrom"]
        s_s, s_e = int(r["start"]), int(r["end"])
        glist = gene_by_chrom.get(chrom, [])
        if not glist:
            continue
        starts = [g[0] for g in glist]
        left = bisect_right(starts, s_s) - 1
        right = bisect_left(starts, s_e)
        best = None
        for idx in range(max(0, left), min(len(glist), right + 2)):
            g_start, g_end, gid = glist[idx]
            if s_s < g_end and s_e > g_start:
                ov = min(s_e, g_end) - max(s_s, g_start)
                if best is None or ov > best[1]:
                    best = (gid, ov)
        if best:
            mapping[r["seq_id"]] = best
    return mapping


def family_of(r) -> str:
    for col in FAM_COLS:
        if r.get(col) == "True":
            return col.replace("is_", "").upper()
    return "?"


def load_fdr() -> dict:
    """gene -> padj"""
    fdr = {}
    with open(FDR_TABLE) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gene = r["gene"].split(".")[0]
            padj = r.get("padj", "")
            if padj not in ("", "NA", "nan"):
                fdr[gene] = float(padj)
    return fdr


def main():
    print("=== 表达/DE 分析 (06_expression) ===")
    core = load_core_set()
    print(f"核心集 : {len(core):,} sORF")

    rows = load_ssp_rows()
    print(f"SSP candidates (Yu11 is_ssp=True): {len(rows):,}")

    genes = load_genes(GFF)
    n_genes = sum(len(v) for v in genes.values())
    print(f"GFF 基因: {n_genes:,}")

    # sORF -> gene
    mapping = map_sorf_to_gene(rows, genes)
    print(f"映射到基因的 SSP sORF: {len(mapping):,} ({len(mapping)/len(rows)*100:.1f}%)")

    # 家族 + 核心集
    fam_of = {r["seq_id"]: family_of(r) for r in rows}
    core_of = {r["seq_id"]: (r["seq_id"] in core) for r in rows}

    # 写 sorf_gene_map.tsv
    with open(OUT / "sorf_gene_map.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sorf_id", "chrom", "start", "end", "strand",
                    "gene", "overlap_bp", "family", "in_core_set",
                    "is_ida_like"])
        for r in rows:
            gid, ov = mapping.get(r["seq_id"], ("", 0))
            w.writerow([r["seq_id"], r["chrom"], r["start"], r["end"],
                        r["strand"], gid, ov, fam_of[r["seq_id"]],
                        core_of[r["seq_id"]], r.get("is_ida_like", "")])
    print(f"✅ sorf_gene_map.tsv  ({len(rows):,} 行)")

    # FPKM
    fpkm = pd.read_csv(FPKM_CSV, index_col=0)
    print(f"FPKM 矩阵: {fpkm.shape[0]} 基因 × {fpkm.shape[1]} 样本")

    yu11_t48_sm = [c for c in fpkm.columns
                   if "Yu11" in c and "4C-" in c and "48h" in c and "Sm" in c]
    yu11_c48_sm = [c for c in fpkm.columns
                   if "Yu11" in c and "CK-" in c and "48h" in c and "Sm" in c]
    base_cols = [c for c in fpkm.columns
                 if "Yu11" in c and "CK-0" in c and "Sm" in c]
    print(f"  4C 48h Sm: {yu11_t48_sm}")
    print(f"  CK 48h Sm: {yu11_c48_sm}")
    print(f"  baseline 0h: {base_cols}")

    # 基因级聚合: gene -> sids
    gene_to_sids = defaultdict(list)
    for r in rows:
        gid, _ = mapping.get(r["seq_id"], ("", 0))
        if gid:
            gene_to_sids[gid].append(r["seq_id"])

    fdr = load_fdr()
    print(f"FDR table (gene-level, revised): {len(fdr):,} genes")

    results = []
    expressed_sorfs = []
    for gid, sids in gene_to_sids.items():
        gid_plus = gid if gid in fpkm.index else (
            f"{gid}.1" if f"{gid}.1" in fpkm.index else None)
        if gid_plus is None:
            continue
        t48 = float(fpkm.loc[gid_plus, yu11_t48_sm].mean())
        c48 = float(fpkm.loc[gid_plus, yu11_c48_sm].mean())
        baseline = float(fpkm.loc[gid_plus, base_cols].mean()) if base_cols else 0.0
        fc = (t48 + 0.1) / (c48 + 0.1)
        log2fc = float(np.log2(fc))
        # expressed sORF 记录
        for sid in sids:
            expressed_sorfs.append({
                "sorf_id": sid, "gene": gid,
                "family": fam_of[sid], "in_core_set": core_of[sid],
                "is_ida_like": next(r.get("is_ida_like", "") for r in rows if r["seq_id"] == sid),
                "t48_stem_fpkm": round(t48, 3), "c48_stem_fpkm": round(c48, 3),
            })
        if max(t48, c48) > 0.5 and abs(log2fc) > 0.5:
            fams = ", ".join(sorted(set(fam_of[s] for s in sids)))
            n_core = sum(1 for s in sids if core_of[s])
            results.append({
                "gene": gid,
                "ssp_ids": "; ".join(sids[:3]),
                "n_sorfs": len(sids),
                "n_core": n_core,
                "ssp_families": fams,
                "log2fc_48h_stem": round(log2fc, 2),
                "t48_stem_fpkm": round(t48, 1),
                "c48_stem_fpkm": round(c48, 1),
                "baseline_0h_fpkm": round(baseline, 1),
                "direction": "UP" if log2fc > 0 else "DOWN",
                "padj": round(fdr.get(gid, float("nan")), 4) if gid in fdr else "NA",
                "fdr_lt_0.05": "Y" if gid in fdr and fdr[gid] < 0.05 else "N",
            })

    results.sort(key=lambda x: abs(x["log2fc_48h_stem"]), reverse=True)
    loose = results
    strict = [r for r in results if abs(r["log2fc_48h_stem"]) > 1.0]
    print(f"\n=== DE 候选 (Yu11 48h Sm, 4C vs CK) ===")
    print(f"宽松 (|log2FC|>0.5): {len(loose)} 基因")
    print(f"严格 (|log2FC|>1.0): {len(strict)} 基因")

    def write_de(fname, rows):
        with open(OUT / fname, "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write_de("de_candidates_loose.tsv", loose)
    write_de("de_candidates_strict.tsv", strict)
    print(f"✅ de_candidates_loose.tsv ({len(loose)}) / strict ({len(strict)})")

    # expressed_sorfs.tsv
    with open(OUT / "expressed_sorfs.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=list(expressed_sorfs[0].keys()))
        w.writeheader()
        w.writerows(expressed_sorfs)
    print(f"✅ expressed_sorfs.tsv ({len(expressed_sorfs):,} sORF)")

    # de_classification.tsv (基因级)
    cls_rows = []
    for r in results:
        cls_rows.append({
            "gene": r["gene"], "log2fc": r["log2fc_48h_stem"],
            "direction": r["direction"], "padj": r["padj"],
            "fdr_lt_0.05": r["fdr_lt_0.05"],
            "n_sorfs": r["n_sorfs"], "n_core": r["n_core"],
            "ssp_families": r["ssp_families"],
            "class": "ssp_de",
        })
    with open(OUT / "de_classification.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(cls_rows[0].keys()))
        w.writeheader()
        w.writerows(cls_rows)
    print(f"✅ de_classification.tsv ({len(cls_rows)})")

    # 汇总统计
    up = sum(1 for r in strict if r["direction"] == "UP")
    down = sum(1 for r in strict if r["direction"] == "DOWN")
    fdr_pass = sum(1 for r in strict if r["fdr_lt_0.05"] == "Y")
    core_hit = sum(1 for r in strict if r["n_core"] > 0)
    fam_dist = defaultdict(int)
    for r in strict:
        for fam in r["ssp_families"].split(", "):
            if fam and fam != "?":
                fam_dist[fam] += 1
    print(f"严格: UP={up} DOWN={down} FDR<0.05={fdr_pass} 含核心集={core_hit}")
    print(f"家族分布: {dict(sorted(fam_dist.items(), key=lambda x: -x[1]))}")

    # Top 候选
    print("\nTop 10 (|log2FC| 排序):")
    for r in strict[:10]:
        print(f"  {r['gene']}: {r['ssp_families']} log2FC={r['log2fc_48h_stem']} "
              f"FDR={r['padj']} core={r['n_core']} {r['direction']}")

    # evidence chain
    n_expr = len(expressed_sorfs)
    n_expr_uniq_gene = len(set(e["gene"] for e in expressed_sorfs))
    chain = f"""# 证据链 (06_expression, Yu11)

Method: identical to the established screen_resistance_ssp.py (stem 48h, 3v3, FPKM means,
log2FC = log2((t+0.1)/(c+0.1)), 严格 |log2FC|>1 & max(t48,c48)>0.5, 宽松 |log2FC|>0.5)。

| 步骤 | 数字 | 说明 |
|:-----|:-----|:-----|
| 1. Layer1 six-frame sORF | 4,451,664 | Layer1 (motif-independent) |
| 2. Layer3 SSP 候选 | {len(rows):,} | motif library is_ssp=True |
| 3. sORF→基因重叠 | {len(mapping):,} ({len(mapping)/len(rows)*100:.1f}%) | GFF gene 坐标重叠, 最大 overlap |
| 4. 有 FPKM 数据 | {n_expr:,} sORF / {n_expr_uniq_gene:,} 基因 | 匹配 FPKM 矩阵 |
| 5. DE 候选 (严格) | {len(strict)} | UP {up} / DOWN {down}, FDR<0.05 {fdr_pass}, 含核心集 {core_hit} |
| 6. DE 候选 (宽松) | {len(loose)} | 同上 0.5 阈值 |

FDR note: gene-level padj reuses the revised deseq2_approx_yu11_sm48.tsv
(count 矩阵基因级, 与 sORF 版本无关; R DESeq2 不可用时的负二项近似 + BH)。

核心集交集: DE 严格候选含 核心集成员的基因 = {core_hit} 个。
"""
    with open(OUT / "evidence_chain.md", "w") as f:
        f.write(chain)
    print(f"\n✅ evidence_chain.md")
    print("\n=== 自检 ===")
    assert len(rows) == 5289, f"SSP 数 != 5289: {len(rows)}"
    assert len(strict) == len([r for r in results if abs(r["log2fc_48h_stem"]) > 1.0])
    assert os.path.exists(OUT / "de_candidates_strict.tsv")
    print("自检通过 ✅")


if __name__ == "__main__":
    main()
