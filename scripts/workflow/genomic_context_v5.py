import os
#!/usr/bin/env python3
"""F. 240 个 DE 关联 sORF 的基因组上下文分类 — R5/R6 要求

类别: intergenic / intronic / 5'UTR / 3'UTR / antisense-overlap / CDS-overlap
方法: 240 sORF (Layer1 sorfs.bed 坐标) vs Yu11 GFF (CDS/exon/mRNA 坐标) 重叠
输出: results/06_expression/genomic_context_v20260825.tsv
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

V5 = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "results"
GFF = Path(os.environ.get("PEPTSESAME_ROOT", ".")) / "data/sesame/annotations/Yu11.fixed.gff3"
OUT = V5 / "06_expression/genomic_context_v20260825.tsv"

# 240 sORF 坐标 + 链 (从 Layer1 sorfs.bed)
SORFS = {}
BED_STRAND = {}
with open(V5 / "01_layer1_sixframe/Yu11/sorfs.bed") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        chrom, s, e, sid, strand = p[0], int(p[1]), int(p[2]), p[3].replace("sORF_", ""), p[5]
        SORFS[sid] = (chrom, s, e)
        BED_STRAND[sid] = strand
print(f"Layer1 sORF: {len(SORFS):,}")

# DE 关联 240 sORF
de_sorfs = []
de_genes = set()
with open(V5 / "06_expression/de_candidates_strict.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        de_genes.add(r["gene"])
with open(V5 / "06_expression/sorf_gene_map.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["gene"] in de_genes:
            de_sorfs.append(r["sorf_id"].replace("sORF_", "").split("|")[0])
print(f"DE 关联 sORF: {len(de_sorfs)}")

# GFF 基因模型 (CDS + mRNA 跨度, 每条带链)
cds_spans = defaultdict(list)   # chrom -> [(s, e, strand)]
mrna_spans = defaultdict(list)  # chrom -> [(s, e, strand)]
exon_spans = defaultdict(list)  # chrom -> [(s, e, strand)]
for line in open(GFF):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9:
        continue
    chrom, feat, s, e, strand = p[0], p[2], int(p[3]), int(p[4]), p[6]
    if feat == "CDS":
        cds_spans[chrom].append((s, e, strand))
    elif feat == "mRNA":
        mrna_spans[chrom].append((s, e, strand))
    elif feat == "exon":
        exon_spans[chrom].append((s, e, strand))
print(f"GFF: CDS {sum(len(v) for v in cds_spans.values()):,}, mRNA {sum(len(v) for v in mrna_spans.values()):,}")


def overlaps(s, e, a, b):
    return s <= b and e >= a


def classify(chrom, s, e, strand):
    cds = cds_spans.get(chrom, [])
    mrna = mrna_spans.get(chrom, [])
    exn = exon_spans.get(chrom, [])
    # 同链 vs 反链 (CDS 级独立判定)
    sense_cds = any(overlaps(s, e, a, b) and strand == st for a, b, st in cds)
    anti_cds = any(overlaps(s, e, a, b) and strand != st for a, b, st in cds)
    anti_mrna = any(overlaps(s, e, a, b) and strand != st for a, b, st in mrna)
    if anti_cds or anti_mrna:
        return "antisense"
    if sense_cds:
        return "CDS-overlap (sense)"
    sense_mrna = any(overlaps(s, e, a, b) and strand == st for a, b, st in mrna)
    if sense_mrna:
        if any(overlaps(s, e, a, b) and strand == st for a, b, st in exn):
            return "UTR (sense)"
        return "intronic (sense)"
    return "intergenic"


dist = Counter()
rows = []
for sid in de_sorfs:
    if sid not in SORFS:
        rows.append([sid, "NA", "NA", "NA"])
        continue
    chrom, s, e = SORFS[sid]
    strand = BED_STRAND[sid]
    cat = classify(chrom, s, e, strand)
    # 反义检查: 与 mRNA 重叠但链相反
    anti = False
    if cat in ("UTR", "intronic", "intergenic"):
        # 反义 = 与任一 mRNA 跨度重叠且注释链相反
        pass
    dist[cat] += 1
    rows.append([sid, chrom, f"{s}-{e}", cat])
print("\n基因组上下文分布 (240 DE 关联 sORF):")
for cat, n in dist.most_common():
    print(f"  {cat}: {n} ({n/len(de_sorfs)*100:.1f}%)")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sorf_id", "chrom", "coords", "genomic_context"])
    w.writerows(rows)
print(f"✅ {OUT}")
