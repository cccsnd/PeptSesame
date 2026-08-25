#!/usr/bin/env python3
"""generate_figures.py — Paper1 图表 数据版 (Figure 2/3/4/5/6/7)

reworked from the earlier version with updated data sources:
- 数据全部来自 results/ (Table2 / TableS1 / 02-03 S3651 / 06 表达 / 07 PTM)
- fig2B classification: the earlier miPEP/AMP columns no longer exist → SSP candidates vs non-SSP
- fig6: 273 DE → 221 DE; summary 文字全部更新为 数字
- fig7: PTM 指纹用 DE 候选 PTM (07 产物)
- 输出 PDF+PNG 到 results/10_figures/
"""
import csv
import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUTDIR = R2 / "10_figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

BLUE = "#2C7FB8"
RED = "#D95F0E"
GREEN = "#31A354"
GRAY = "#999999"
PURPLE = "#756BB1"
FAMILIES = ["CLE", "RALF", "CEP", "PSK", "PSY1", "IDA", "EPFL", "RGF"]

sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 11, "figure.dpi": 150})


def save(fig, name):
    fig.savefig(f"{OUTDIR}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ {name} (pdf+png)")


def load_table2():
    """合并家族计数 + 基因组统计"""
    fam = pd.read_csv(R2 / "09_tables/table2_family_counts.csv")
    stats = pd.read_csv(R2 / "09_tables/TableS1_18genome_stats.tsv", sep="\t")
    stats = stats.rename(columns={"Species": "species"})
    df = fam.merge(stats, on="species", suffixes=("", "_y"))
    df = df.rename(columns={
        "Size_Mb": "genome_size_mb",
        "Density_per_Mb": "sorf_density",
        "SSPs": "n_ssp",
    })
    return df


def load_confidence(path):
    counts = Counter()
    with open(path) as f:
        header = f.readline().strip().split("\t")
        idx = header.index("confidence")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > idx:
                counts[parts[idx].strip()] += 1
    return counts


def load_length_dist(bed_path):
    l10 = l20 = l50 = 0
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            aa = int(parts[8])
            if 10 <= aa < 20:
                l10 += 1
            elif 20 <= aa < 50:
                l20 += 1
            elif 50 <= aa <= 100:
                l50 += 1
    return l10, l20, l50


def sp_color(s):
    if s in ("Yu11", "S3651"):
        return RED
    if s in ("14G01", "14G02", "K16", "ken1", "ken8"):
        return GREEN
    return BLUE


# ═══════════════ Figure 2 ═══════════════
def fig2(table2):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    ax = axes[0, 0]
    ax.set_title("A Scoring Confidence (S3651)", fontweight="bold")
    conf = load_confidence(R2 / "02_layer2_scoring/S3651/scored_sorfs.tsv")
    sizes = [conf.get("high", 0), conf.get("medium", 0), conf.get("low", 0)]
    ax.pie(sizes, labels=["HIGH", "MEDIUM", "LOW"],
           autopct="%1.1f%%", colors=[GREEN, RED, GRAY], startangle=90)
    ax.text(0, 0, f"{sum(sizes):,}\ntotal", ha="center", va="center",
            fontsize=10, fontweight="bold")

    ax = axes[0, 1]
    ax.set_title("B sORF Classification (S3651)", fontweight="bold")
    n_ssp = 0
    n_total = 0
    with open(R2 / "03_layer3_classify/S3651/classified_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            n_total += 1
            if r["is_ssp"] == "True":
                n_ssp += 1
    cat_data = {"SSP candidates": n_ssp, "Non-SSP": n_total - n_ssp}
    bars = ax.barh(list(cat_data.keys()), list(cat_data.values()),
                   color=[RED, GRAY])
    for bar, val in zip(bars, cat_data.values()):
        ax.text(bar.get_width() + 20000, bar.get_y() + bar.get_height()/2,
                f"{val:,} ({val/n_total*100:.2f}%)", va="center", fontsize=8)
    ax.set_xlabel("Count")

    ax = axes[0, 2]
    ax.set_title("C Peptide Length (S3651)", fontweight="bold")
    l10, l20, l50 = load_length_dist(R2 / "01_layer1_sixframe/S3651/sorfs.bed")
    len_data = {"10-20aa": l10, "20-50aa": l20, "50-100aa": l50}
    ax.pie(len_data.values(), labels=len_data.keys(),
           autopct=lambda p: f"{p:.1f}%\n({int(p*sum(len_data.values())/100):,})",
           colors=[RED, BLUE, GREEN], startangle=90, textprops={"fontsize": 8})

    ax = axes[1, 0]
    ax.set_title("D SSP Family Members", fontweight="bold")
    s3651_row = table2[table2["species"] == "S3651"].iloc[0]
    yu11_row = table2[table2["species"] == "Yu11"].iloc[0]
    x = np.arange(len(FAMILIES))
    w = 0.38
    ax.bar(x - w/2, [s3651_row[f] for f in FAMILIES], w, label="S3651 (wild)", color=RED, alpha=0.85)
    ax.bar(x + w/2, [yu11_row[f] for f in FAMILIES], w, label="Yu11 (cult.)", color=BLUE, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(FAMILIES, fontsize=9)
    ax.set_ylabel("Count")
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ax.set_title("E Wild vs Cultivated Sesame", fontweight="bold")
    metrics = ["Genome\nSize", "sORF\nDensity", "SSP\nCount"]
    s3651_v = [s3651_row["genome_size_mb"], s3651_row["sorf_density"], s3651_row["n_ssp"]]
    yu11_v = [yu11_row["genome_size_mb"], yu11_row["sorf_density"], yu11_row["n_ssp"]]
    x = np.arange(len(metrics))
    w = 0.3
    ax.bar(x - w/2, s3651_v, w, label="S3651 (wild)", color=RED, alpha=0.8)
    ax.bar(x + w/2, yu11_v, w, label="Yu11 (cult.)", color=BLUE, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.legend(fontsize=9)

    axes[1, 2].axis("off")  # cross-species density panel moved to Figure 4A

    plt.tight_layout()
    save(fig, "Figure2_sorf_catalog")


# ═══════════════ Figure 3 ═══════════════
def fig3(table2):
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))

    mat = table2.copy()
    for f in FAMILIES:
        mat[f] = mat[f] / mat["genome_size_mb"] * 100
    mat = mat.set_index("species")[FAMILIES]

    ax = axes[0]
    sns.heatmap(mat, ax=ax, cmap="YlOrRd", annot=True, fmt=".1f",
                linewidths=0.5, cbar_kws={"label": "SSPs per 100 Mb"},
                annot_kws={"fontsize": 7})
    ax.set_title("A SSP Family Density (18 genomes, per 100 Mb)", fontweight="bold")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    ax = axes[1]
    totals = table2.sort_values("n_ssp")
    bars = ax.barh(totals["species"], totals["n_ssp"],
                   color=[sp_color(s) for s in totals["species"]])
    for bar, val in zip(bars, totals["n_ssp"]):
        ax.text(bar.get_width() + max(totals["n_ssp"])*0.01,
                bar.get_y() + bar.get_height()/2, f"{val:,}", va="center", fontsize=8)
    ax.set_xlabel("SSP Candidates")
    ax.set_title("B SSP Count by Species", fontweight="bold")

    plt.tight_layout()
    save(fig, "Figure3_ssp_families")


# ═══════════════ Figure 4 ═══════════════
def fig4(table2):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    df = table2
    colors = [sp_color(s) for s in df["species"]]

    ax = axes[0]
    ax.scatter(df["genome_size_mb"], df["sorf_density"], c=colors, s=110,
               alpha=0.85, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(row["species"], (row["genome_size_mb"], row["sorf_density"]),
                    fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("Genome Size (Mb)")
    ax.set_ylabel("sORF Density (/Mb)")
    ax.set_title("A sORF Density vs Genome Size (18 genomes)", fontweight="bold")

    ax = axes[1]
    ax.scatter(df["genome_size_mb"], df["n_ssp"], c=colors, s=110, alpha=0.85, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(row["species"], (row["genome_size_mb"], row["n_ssp"]),
                    fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("Genome Size (Mb)")
    ax.set_ylabel("SSP Candidates")
    ax.set_title("B SSP Candidates vs Genome Size", fontweight="bold")

    plt.tight_layout()
    save(fig, "Figure4_cross_species")


# ═══════════════ Figure 5 ═══════════════
def fig5():
    """DE Top 基因表达热图 (Yu11 茎 48h 4C vs CK + baseline)"""
    de = pd.read_csv(R2 / "06_expression/de_candidates_strict.tsv", sep="\t")
    fpkm = pd.read_csv(R2 / "00_inputs/rnaseq/FPKM_matrix.csv", index_col=0)

    top = de.reindex(de["log2fc_48h_stem"].abs().sort_values(ascending=False).index).head(20)
    cols = ["Yu11-4C-48h-Sm1", "Yu11-4C-48h-Sm2", "Yu11-4C-48h-Sm3",
            "Yu11-CK-48h-Sm1", "Yu11-CK-48h-Sm2", "Yu11-CK-48h-Sm3"]
    rows = []
    for g in top["gene"]:
        gid = g if g in fpkm.index else f"{g}.1"
        if gid in fpkm.index:
            rows.append(fpkm.loc[gid, cols])
    expr = pd.DataFrame(rows, index=top["gene"].iloc[:len(rows)])

    fig, ax = plt.subplots(figsize=(10, max(6, len(expr)*0.35)))
    sns.heatmap(np.log2(expr + 1), ax=ax, cmap="RdYlBu_r", center=0,
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "log2(FPKM+1)"})
    ax.set_title("Top 20 DE SSP-candidate genes (Yu11 stem, 48h, 4C vs CK)",
                 fontweight="bold")
    plt.tight_layout()
    save(fig, "Figure5_expression_heatmap")


# ═══════════════ Figure 6 ═══════════════
def fig6():
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    ax = axes[0]
    ax.axis("off")
    cands = pd.read_csv(R2 / "06_expression/de_candidates_strict.tsv", sep="\t")
    top = cands.sort_values("log2fc_48h_stem", ascending=False).head(10)
    row_data = [[r["gene"], r["ssp_families"], f"{r['log2fc_48h_stem']:.2f}",
                 f"{r['t48_stem_fpkm']:.0f}", f"{r['c48_stem_fpkm']:.0f}",
                 r["direction"]] for _, r in top.iterrows()]
    table = ax.table(cellText=row_data,
                     colLabels=["Gene", "Family", "log2FC", "T48Sm", "C48Sm", "Dir"],
                     loc="center", cellLoc="center",
                     colWidths=[0.16, 0.12, 0.14, 0.14, 0.14, 0.1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)
    ax.set_title("A Top 10 Induced SSP Candidates (Stem, 48h)", fontweight="bold", pad=20)

    ax = axes[1]
    ax.axis("off")
    table2 = load_table2()
    dens = table2["sorf_density"]
    cv = dens.std() / dens.mean() * 100
    stats_text = (
        "PeptSesame Output Summary\n"
        "══════════════════════════════════════════\n\n"
        "Layer 1: Six-frame translation\n"
        "  Yu11:   4,451,664 sORFs (CDS-filtered)\n"
        "  S3651:  8,581,265 sORFs\n\n"
        "Layer 2: Scoring (rule-based, 6 channels)\n"
        "  weight sensitivity: ρ > 0.99\n\n"
        "Layer 3: Classification (motifs)\n"
        "  SSP candidates: 5,289 (Yu11)\n"
        "                   4,862 (S3651)\n\n"
        "Layer 4: Function & Expression \n"
        "  1,093 SSPs overlap genes (20.7%)\n"
        "  942 have RNA-seq expression (17.8%)\n"
        "  221 DE genes in stem at 48h\n"
        "    UP: 106, DOWN: 115 (30 FDR<0.05)\n\n"
        "Cross-species (18 genomes):\n"
        f"  sORF density: 12,700-15,800/Mb (CV {cv:.1f}%)\n"
        "  All 8 SSP families detected\n"
        "  IDA DE candidates removed (motif fix)"
    )
    ax.text(0.03, 0.95, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment="top",
            fontfamily="monospace", linespacing=1.5)
    ax.set_title("B Pipeline Summary ", fontweight="bold")

    plt.tight_layout()
    save(fig, "Figure6_functional_candidates")


# ═══════════════ Figure 7 ═══════════════
def fig7():
    """SSP 家族 PTM 指纹 (DE 候选 PTM, per-member 位点密度)"""
    ptm = pd.read_csv(R2 / "07_ptm/plantptm_de_candidates.tsv", sep="\t")
    # member counts: family sORF counts (Yu11)
    n_member = Counter()
    with open(R2 / "03_layer3_classify/Yu11/classified_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] == "True":
                for fam in r["ssp_families"].split(";"):
                    if fam:
                        n_member[fam] += 1

    PTM_TYPES = sorted(ptm["ptm_type"].unique())
    PTM_LABEL = {"pho": "Phospho", "Kac": "K-acetyl", "Kcr": "K-crotonyl",
                 "Khib": "K-2-hydroxyisobutyryl", "Kmal": "K-malonyl",
                 "Ksucc": "K-succinyl", "Kub": "K-ubiquitin",
                 "Ngly": "N-glycosyl", "Sacy": "S-acyl"}
    fams = sorted(set(ptm["family"]))
    mat = pd.DataFrame(0.0, index=fams, columns=PTM_TYPES)
    for _, r in ptm.iterrows():
        if r["family"] in mat.index and r["ptm_type"] in mat.columns:
            mat.loc[r["family"], r["ptm_type"]] += 1
    mat = mat / mat.sum(axis=1).replace(0, np.nan).values[:, None] * 100
    mat = mat.fillna(0)

    fig, ax = plt.subplots(figsize=(11, 5))
    sns.heatmap(mat, ax=ax, cmap="YlOrRd", annot=True, fmt=".0f",
                linewidths=0.5, cbar_kws={"label": "% of PTM sites"},
                annot_kws={"fontsize": 8})
    ax.set_yticklabels([f"{f} (n={n_member[f]:,})" for f in fams],
                       rotation=0, fontsize=9)
    ax.set_title("SSP Family PTM Composition (DE candidates, PlantPTM high-confidence)",
                 fontweight="bold")
    plt.tight_layout()
    save(fig, "Figure7_ptm_fingerprint")


if __name__ == "__main__":
    t2 = load_table2()
    print(f"Table2: {len(t2)} species × {len(t2.columns)} columns")
    fig2(t2)
    fig3(t2)
    fig4(t2)
    fig5()
    fig6()
    fig7()
    print(f"\nAll figures (pdf+png) saved to {OUTDIR}/")
