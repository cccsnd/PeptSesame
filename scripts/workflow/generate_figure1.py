#!/usr/bin/env python3
"""generate_figure1.py — Figure 1: pipeline schematic + real-data statistics


- Panel A: 专业管线流程图 (matplotlib FancyBboxPatch + 箭头)
- Panel B: 真实基因组统计 (修正 Yu11: 4,451,664 sORF / 4,456 SSP)
- Panel C: 真实sORF长度分布 (从BED计算, 非随机!)
- Panel D: 真实染色体分布 (从BED计算)
"""
import csv, sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import os as _os
from pathlib import Path as _Path
ROOT = _Path(_os.environ.get("PEPTSESAME_ROOT", "."))
OUTDIR = ROOT / "results/10_figures"
OUTDIR.mkdir(parents=True, exist_ok=True)
BLUE = "#2C7FB8"; RED = "#D95F0E"; GREEN = "#31A354"; GRAY = "#999999"; LIGHT = "#E8F0F8"; GOLD = "#B8860B"

# ─── 真实数据加载 ───────────────────────────────────────────
def load_bed_stats(bed_path):
    lengths = []
    chr_counts = Counter()
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            nt = int(parts[7])
            if 30 <= nt <= 300:
                lengths.append(nt)
            chrom = parts[0].split("#")[-1].replace("Chr", "chr")
            chr_counts[chrom] += 1
    return lengths, chr_counts

s3651_bed = ROOT / "results/01_layer1_sixframe/S3651/sorfs.bed"
yu11_bed = ROOT / "results/01_layer1_sixframe/Yu11/sorfs.bed"
s3651_len, s3651_chr = load_bed_stats(s3651_bed)
yu11_len, yu11_chr = load_bed_stats(yu11_bed)
print(f"S3651: {len(s3651_len):,} sORFs, {len(s3651_chr)} chroms")
print(f"Yu11:  {len(yu11_len):,} sORFs, {len(yu11_chr)} chroms")

# ─── Figure ─────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], width_ratios=[1, 1],
                      hspace=0.35, wspace=0.3)

# ═══ Panel A: 管线流程图 ═══
axA = fig.add_subplot(gs[0, 0])
axA.set_xlim(0, 10); axA.set_ylim(0, 10)
axA.axis("off")
axA.set_title("A  PeptSesame Pipeline", loc="left", fontweight="bold", fontsize=13)

def box(ax, x, y, w, h, text, fc=LIGHT, ec=BLUE, fs=9, bold=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                       fc=fc, ec=ec, lw=1.5)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fs, fontweight="bold" if bold else "normal")

def arrow(ax, x1, y1, x2, y2, color=BLUE):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                 arrowstyle="-|>", mutation_scale=18, color=color, lw=1.8))

# Input
box(axA, 3.2, 9.1, 3.6, 0.6, "Genome FASTA (+GFF, +RNA-seq)", fc="#FFF3E0", ec=RED, fs=9, bold=True)
arrow(axA, 5.0, 9.1, 5.0, 8.6)

# Layer 1
box(axA, 2.2, 7.4, 5.6, 1.2,
    "Layer 1: Six-frame translation\nORF extraction (30–300 nt) + CDS filtering",
    fc=LIGHT, ec=BLUE, fs=9)
arrow(axA, 5.0, 7.4, 5.0, 6.9)

# Layer 2
box(axA, 2.2, 5.7, 5.6, 1.2,
    "Layer 2: Multi-evidence scoring\n2 core channels + orthogonal evidence layers",
    fc=LIGHT, ec=BLUE, fs=9)
arrow(axA, 5.0, 5.7, 5.0, 5.2)

# Layer 3
box(axA, 2.2, 4.0, 5.6, 1.2,
    "Layer 3: Classification\n8 SSP families + IDA-like (motifs)",
    fc=LIGHT, ec=BLUE, fs=9)
arrow(axA, 5.0, 4.0, 5.0, 3.5)

# Layer 4
box(axA, 2.2, 2.3, 5.6, 1.2,
    "Layer 4: Functional annotation\nGO / KEGG / stress association / receptor pairing",
    fc=LIGHT, ec=BLUE, fs=9)
arrow(axA, 5.0, 2.3, 5.0, 1.8)

# Output
box(axA, 2.8, 0.9, 4.4, 0.7,
    "sORF catalog + scoring + classification", fc="#E8F5E9", ec=GREEN, fs=9, bold=True)

# Right-side: 18 genomes badge
box(axA, 8.3, 4.0, 1.5, 2.4, "18\nplant\ngenomes", fc="#F3E5F5", ec="#7B1FA2", fs=10, bold=True)
arrow(axA, 7.8, 5.2, 8.3, 5.2, color="#7B1FA2")

# Right-side (added): 验证短名单 + 外部证据交叉标注
box(axA, 7.9, 2.2, 1.9, 1.4,
    "Validation\nshortlists\nA 3 / B 55 / C 240", fc="#FFF8E1", ec=GOLD, fs=8, bold=True)
arrow(axA, 7.2, 1.4, 7.9, 2.6, color=GOLD)
box(axA, 7.9, 0.4, 1.9, 1.4,
    "External\nevidence\ncross-referencing", fc="#E0F2F1", ec=GREEN, fs=7.5, bold=True)
arrow(axA, 7.2, 1.2, 7.9, 1.0, color=GREEN)

# ═══ Panel B: 基因组统计 ═══
axB = fig.add_subplot(gs[0, 1])
axB.axis("off")
axB.set_title("B  Sesame T2T Genome Statistics", loc="left", fontweight="bold", fontsize=13)
stats_text = (
    "                        S3651 (wild)    Yu11 (cult.)\n"
    "─────────────────────────────────────────────────\n"
    "Species:                S. alatum       S. indicum\n"
    "Genome size:            552 Mb          305 Mb\n"
    "Chromosomes:            13              13\n"
    "sORFs:                  8,581,265       4,451,664\n"
    "sORF density:           15,552/Mb       14,592/Mb\n"
    "SSP candidates:         4,862           5,289\n"
    "RNA-seq samples:        —               132 (M133+Yu11)\n"
    "DE candidates (48h):    —               221 (106 UP / 115 DOWN)"
)
axB.text(0.02, 0.95, stats_text, transform=axB.transAxes,
         fontsize=11, verticalalignment="top", fontfamily="monospace", linespacing=1.7)

# ═══ Panel C: 真实长度分布 ═══
axC = fig.add_subplot(gs[1, 0])
axC.set_title("C  sORF Length Distribution (real data)", fontweight="bold", fontsize=12)
bins = np.arange(30, 305, 10)
axC.hist(s3651_len, bins=bins, color=RED, alpha=0.7, label=f"S3651 (n={len(s3651_len):,})")
axC.hist(yu11_len, bins=bins, color=BLUE, alpha=0.55, label=f"Yu11 (n={len(yu11_len):,})")
mean_s = np.mean(s3651_len); mean_y = np.mean(yu11_len)
axC.axvline(mean_s, color=RED, linestyle="--", lw=1, label=f"S3651 mean {mean_s:.0f} nt")
axC.axvline(mean_y, color=BLUE, linestyle="--", lw=1, label=f"Yu11 mean {mean_y:.0f} nt")
axC.set_xlabel("sORF length (nt)"); axC.set_ylabel("Count")
axC.legend(fontsize=8)
axC.set_yscale("log")

# ═══ Panel D: 染色体分布 ═══
axD = fig.add_subplot(gs[1, 1])
axD.set_title("D  sORF Distribution per Chromosome", fontweight="bold", fontsize=12)
all_chroms = sorted(set(s3651_chr) | set(yu11_chr), key=lambda c: int(c.replace("chr", "")))
x = np.arange(len(all_chroms))
w = 0.38
s3651_vals = [s3651_chr.get(c, 0) / 1e6 for c in all_chroms]
yu11_vals = [yu11_chr.get(c, 0) / 1e6 for c in all_chroms]
axD.bar(x - w/2, s3651_vals, w, label="S3651", color=RED, alpha=0.8)
axD.bar(x + w/2, yu11_vals, w, label="Yu11", color=BLUE, alpha=0.8)
axD.set_xticks(x); axD.set_xticklabels(all_chroms, fontsize=7, rotation=45)
axD.set_xlabel("Chromosome"); axD.set_ylabel("sORFs (millions)")
axD.legend(fontsize=9)

plt.savefig(f"{OUTDIR}/Figure1_pipeline_overview.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{OUTDIR}/Figure1_pipeline_overview.pdf", bbox_inches="tight")
plt.close()
print(f"✅ Saved Figure1_pipeline_overview (png+pdf)")
print(f"   S3651 mean length: {mean_s:.1f} nt, Yu11: {mean_y:.1f} nt")
