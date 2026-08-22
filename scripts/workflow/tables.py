import os
#!/usr/bin/env python3
"""tables.py — 09_tables 核心表生成 (current pipeline)

生成:
  - table2_family_counts.csv   Table 2: 18 物种 SSP 家族计数 (源: family_counts.tsv)
  - TableS1_18genome_stats.tsv 18 物种基因组统计 (Size/sORF 数/密度/SSP/家族)
  - TableS5_DE_candidates.tsv  DE 候选 (严格, 221)
输出到 results_v2/09_tables/
"""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results_v2"
OUT = R2 / "09_tables"
OUT.mkdir(parents=True, exist_ok=True)

GENOMES = R2 / "00_inputs/genomes"
FAM_COUNTS = OUT / "family_counts.tsv"
DE_STRICT = R2 / "06_expression/de_candidates_strict.tsv"

FAM_ORDER = ["CLE", "RALF", "CEP", "PSK", "PSY1", "IDA", "EPFL", "RGF"]


def genome_size_mb(fasta: Path) -> float:
    """FASTA 总碱基数 (Mb)"""
    n = 0
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                continue
            n += len(line.strip())
    return n / 1e6


def count_sorfs(species: str) -> int:
    bed = R2 / f"01_layer1_sixframe/{species}/sorfs.bed"
    if not bed.exists():
        return -1
    return sum(1 for _ in open(bed)) - 1


def main():
    # 读家族计数
    fam_rows = {}
    with open(FAM_COUNTS) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fam_rows[r["species"]] = r

    # species -> fasta 名映射
    fasta_files = {p.stem: p for p in GENOMES.glob("*.fasta")}

    stats = []
    for sp, r in fam_rows.items():
        fa = fasta_files.get(sp) or fasta_files.get(sp.replace("_", "")) \
            or next((p for s, p in fasta_files.items() if s.startswith(sp)), None)
        if fa is None:
            print(f"⚠️ {sp}: 无匹配 fasta, 跳过 Size")
            size = float("nan")
        else:
            size = genome_size_mb(fa)
        n_sorf = count_sorfs(sp)
        density = n_sorf / size if size == size and size > 0 else float("nan")
        stats.append({
            "Species": sp, "Size_Mb": round(size, 2),
            "sORFs": n_sorf, "Density_per_Mb": round(density, 1),
            "SSPs": int(r["ssp_total"]),
            "ida_like": int(r["ida_like"]),
            **{f: int(r[f]) for f in FAM_ORDER},
        })

    # Table 2 (CSV)
    with open(OUT / "table2_family_counts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["species", "ssp_total", "ida_like"] + FAM_ORDER)
        w.writeheader()
        for sp, r in fam_rows.items():
            w.writerow({"species": sp, "ssp_total": r["ssp_total"],
                        "ida_like": r["ida_like"], **{f: r[f] for f in FAM_ORDER}})
    print(f"✅ table2_family_counts.csv")

    # Table S1
    with open(OUT / "TableS1_18genome_stats.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["Species", "Size_Mb", "sORFs",
                                       "Density_per_Mb", "SSPs", "ida_like"] + FAM_ORDER)
        w.writeheader()
        for s in stats:
            w.writerow(s)
    print(f"✅ TableS1_18genome_stats.tsv")
    for s in stats:
        print(f"  {s['Species']:12s} {s['Size_Mb']:8.2f}Mb  {s['sORFs']:>10,} sORF  "
              f"{s['Density_per_Mb']:8.1f}/Mb  SSP={s['SSPs']:,}")

    # Table S5 (DE 候选)
    with open(OUT / "TableS5_DE_candidates.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["Gene", "SSP_Families", "n_sORFs",
                                       "log2FC_48h_Stem", "T48_FPKM", "C48_FPKM",
                                       "Baseline_0h_FPKM", "Direction", "padj",
                                       "FDR_lt_0.05"])
        w.writeheader()
        n = 0
        with open(DE_STRICT) as fsrc:
            for r in csv.DictReader(fsrc, delimiter="\t"):
                w.writerow({
                    "Gene": r["gene"], "SSP_Families": r["ssp_families"],
                    "n_sORFs": r["n_sorfs"],
                    "log2FC_48h_Stem": r["log2fc_48h_stem"],
                    "T48_FPKM": r["t48_stem_fpkm"],
                    "C48_FPKM": r["c48_stem_fpkm"],
                    "Baseline_0h_FPKM": r["baseline_0h_fpkm"],
                    "Direction": r["direction"], "padj": r["padj"],
                    "FDR_lt_0.05": r["fdr_lt_0.05"],
                })
                n += 1
    print(f"✅ TableS5_DE_candidates.tsv ({n} 候选)")


if __name__ == "__main__":
    main()
