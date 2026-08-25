import os
#!/usr/bin/env python3
"""ptm.py — DE 候选 PTM 注释 (07_ptm, 补充级)

Strategy: sORF set of the 221 strict DE genes → intersect with earlier PlantPTM high-confidence results by sORF id
匹配复用 (PlantPTM 位点注释与 sORF 版本无关, Layer1 产物相同, sORF id 一致)。
Earlier results:
  - plantptm_annotation_yu11_273cand_highconf.tsv (273 earlier DE genes)
  - plantptm_annotation_yu11_allssp_highconf.tsv (Yu11 全部 SSP)
sORFs from genes beyond the original 273 that are absent from the earlier table are flagged as missing.

输出: results/07_ptm/
  - plantptm_de_candidates.tsv   DE sORF × PTM 位点 (高置信)
  - ptm_coverage.md              覆盖统计
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results"
OUT = R2 / "07_ptm"
OUT.mkdir(parents=True, exist_ok=True)

V1_273 = ROOT / "results/database/plantptm_annotation_yu11_273cand_highconf.tsv"
V1_ALL = ROOT / "results/database/plantptm_annotation_yu11_allssp_highconf.tsv"
DE_STRICT = R2 / "06_expression/de_candidates_strict.tsv"
SORF_MAP = R2 / "06_expression/sorf_gene_map.tsv"


def load_ptm(path) -> dict:
    """sorf_id -> [rows]"""
    d = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["Protein"]].append(r)
    return d


def main():
    ptm273 = load_ptm(V1_273)
    ptm_all = load_ptm(V1_ALL)
    print(f"earlier PTM sites: 273cand table {sum(len(v) for v in ptm273.values()):,} sites "
          f"({len(ptm273):,} sORF); allssp 表 {sum(len(v) for v in ptm_all.values()):,} "
          f"位点 ({len(ptm_all):,} sORF)")

    # DE 基因 -> sORF
    de_genes = []
    with open(DE_STRICT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            de_genes.append(r["gene"])
    de_set = set(de_genes)

    sorf_of_gene = defaultdict(list)
    with open(SORF_MAP) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["gene"] in de_set:
                sorf_of_gene[r["gene"]].append(
                    (r["sorf_id"], r["family"], r["in_core_set"]))

    all_de_sorfs = [s for g in de_set for s, _, _ in sorf_of_gene[g]]
    print(f"strict DE: {len(de_set)} genes, {len(all_de_sorfs):,} sORFs")

    # 覆盖统计
    hit_273 = sum(1 for s in all_de_sorfs if s in ptm273)
    hit_all = sum(1 for s in all_de_sorfs if s in ptm_all)
    covered = set(all_de_sorfs) & (set(ptm273) | set(ptm_all))
    missing = set(all_de_sorfs) - covered
    print(f"PTM 覆盖: 273表 {hit_273}, allssp表 {hit_all}, 并集覆盖 {len(covered):,} "
          f"({len(covered)/len(all_de_sorfs)*100:.1f}%), 缺失 {len(missing):,}")

    # 输出: DE sORF × PTM (并集, 去重按 sorf+type+pos)
    out_rows = []
    seen = set()
    for s in all_de_sorfs:
        for r in ptm273.get(s, []) + ptm_all.get(s, []):
            key = (s, r["PTM type"], r["Position"])
            if key in seen:
                continue
            seen.add(key)
            # 用 的家族标注
            fam = "?"
            for g in de_set:
                for sid, f, _ in sorf_of_gene[g]:
                    if sid == s:
                        fam = f
                        break
                if fam != "?":
                    break
            out_rows.append({
                "sorf_id": s, "ptm_type": r["PTM type"],
                "position": r["Position"], "residue": r["Residue"],
                "score": r["Score"], "confidence": r["Confidence"],
                "ptm_cn": r["ptm_cn"], "family": fam,
            })

    with open(OUT / "plantptm_de_candidates.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["sorf_id", "ptm_type", "position",
                                       "residue", "score", "confidence",
                                       "ptm_cn", "family"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"✅ plantptm_de_candidates.tsv ({len(out_rows):,} 位点)")

    fam_ptm = Counter()
    for r in out_rows:
        fam_ptm[r["family"]] += 1
    print(f"家族×PTM: {dict(fam_ptm)}")

    # 覆盖报告
    n_missing_genes = len({s.split("_")[1] for s in missing}) if missing else 0
    md = f"""# PTM 覆盖报告 (DE 候选, 补充级)

- 严格 DE 候选: {len(de_set)} 基因 / {len(all_de_sorfs):,} sORF
- PTM 高置信位点覆盖: {len(covered):,} sORF ({len(covered)/len(all_de_sorfs)*100:.1f}%)
- missing: {len(missing):,} sORF (flagged)
- 输出位点: {len(out_rows):,} (去重 sorf×type×pos)

Source: earlier PlantPTM high-confidence results (273cand + allssp tables, same Layer1 products,
sORF id 一致, 位点注释可直接复用)。
"""
    with open(OUT / "ptm_coverage.md", "w") as f:
        f.write(md)
    print(f"✅ ptm_coverage.md")


if __name__ == "__main__":
    main()
