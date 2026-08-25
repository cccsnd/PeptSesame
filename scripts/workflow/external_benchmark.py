#!/usr/bin/env python3
"""external_benchmark.py - independent external benchmark (Wang et al. 2020 peptides)

数据: Wang et al. 2020 Mol Plant 补充表 (documents/maize_ara/)
  - mmc13 = Table S11: Arabidopsis NCPs (MS 验证的非经典肽, 1,862 条)
  - mmc14 = Table S12: Arabidopsis CPs  (MS 验证的常规肽, 2,365 条)
基准 1 (坐标级召回): MS 验证肽的基因组坐标 vs PeptSesame Arabidopsis
  sORF 坐标 (results/03/Arabidopsis) → 被任意 sORF 捕获比例 + 被
  is_ssp=True 候选捕获比例 (nuclear 染色体 1-5; 排除 Mt/Pt)。
基准 2 (序列级 motif 命中): 八家族 motif 在 CP/NCP 序列上的命中率
  (CP 中大部分不是八家族成员, 命中比例是 motif 特异性的独立参考)。
基准 3 (负集 FPR): Yu11 注释蛋白 60-aa 滑窗片段 (非 SSP 蛋白背景,
  组成匹配) 的八家族 motif 命中率 vs 随机 60-aa 肽。

输出: results/08_benchmark/external_benchmark.tsv + 控制台摘要
"""
import csv
import re
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
sys.path.insert(0, str(ROOT))
from pipeline.layer3_classify.motif_profiles import SSP_MOTIFS

FAMS = [f for f in SSP_MOTIFS if f != "IDA_LIKE"]


def read_peptide_table(xlsx, sheet_rows_min=4, cols=(0, 1, 2, 3, 4)):
    """读 xlsx 肽表: (seq, chr, start, end, strand)"""
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, read_only=True)
    ws = wb.worksheets[0]
    recs = []
    for row in ws.iter_rows(min_row=sheet_rows_min, values_only=True):
        if not row or not row[0]:
            continue
        seq = str(row[0]).strip().upper()
        chrom = str(row[1]).strip()
        if chrom in ("Mt", "Pt", "None", ""):
            continue
        try:
            start, end = int(row[2]), int(row[3])
        except (TypeError, ValueError):
            continue
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
            continue
        recs.append((seq, chrom, start, end, str(row[4] or "")))
    return recs


def load_arabidopsis_sorfs():
    """chrom -> [(start, end, is_ssp)]"""
    from collections import defaultdict
    idx = defaultdict(list)
    path = ROOT / "results/03_layer3_classify/Arabidopsis/classified_sorfs.tsv"
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            idx[r["chrom"]].append((int(r["start"]), int(r["end"]),
                                    r["is_ssp"] == "True"))
    for c in idx:
        idx[c].sort()
    return idx


def overlaps(idx, chrom, start, end):
    from bisect import bisect_left
    lst = idx.get(chrom, [])
    if not lst:
        return None
    starts = [x[0] for x in lst]
    i = bisect_left(starts, start)
    best = None
    for j in range(max(0, i - 1), min(len(lst), i + 2)):
        s, e, ssp = lst[j]
        if s < end and e > start:
            ov = min(e, end) - max(s, start)
            if best is None or ov > best[0]:
                best = (ov, ssp)
    return best


def motif_hits(seq):
    return [f for f in FAMS if re.search(SSP_MOTIFS[f], seq)]


def main():
    cp = read_peptide_table(ROOT / "documents/maize_ara/mmc14.xlsx")
    ncp = read_peptide_table(ROOT / "documents/maize_ara/mmc13.xlsx")
    print(f"CP: {len(cp)} 条 (核 1-5), NCP: {len(ncp)} 条")

    idx = load_arabidopsis_sorfs()
    out = []

    for label, recs in (("CP", cp), ("NCP", ncp)):
        n_any = n_ssp = 0
        mh_any = 0
        len_ok = 0
        for seq, chrom, start, end, strand in recs:
            if 10 <= len(seq) <= 100:
                len_ok += 1
            if motif_hits(seq):
                mh_any += 1
            ov = overlaps(idx, chrom, start, end)
            if ov:
                n_any += 1
                if ov[1]:
                    n_ssp += 1
        out.append({
            "set": label, "n": len(recs),
            "len_10_100": len_ok,
            "covered_by_any_sorf": n_any,
            "covered_by_ssp_candidate": n_ssp,
            "motif_hit_rate": mh_any,
        })
        print(f"{label}: {len(recs)} 条 | 长度10-100aa: {len_ok} | "
              f"被任意 sORF 覆盖: {n_any} ({n_any/len(recs)*100:.1f}%) | "
              f"被 SSP 候选覆盖: {n_ssp} ({n_ssp/len(recs)*100:.1f}%) | "
              f"序列 motif 命中: {mh_any} ({mh_any/len(recs)*100:.1f}%)")

    # 负集: Yu11 注释蛋白 60aa 滑窗
    gff = ROOT / "results/00_inputs/annotations/Yu11.gff3"
    fa = ROOT / "results/00_inputs/genomes/Yu11.fasta"
    from collections import defaultdict
    cds = defaultdict(list)
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            cds[p[0]].append((int(p[3]), int(p[4]), p[6]))
    chrom_seqs = {}
    cur = None
    buf = []
    with open(fa) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    chrom_seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur:
            chrom_seqs[cur] = "".join(buf)
    comp = str.maketrans("ACGT", "TGCA")
    CODON = {}
    _bases = "TCAG"
    _aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    for i, b1 in enumerate(_bases):
        for j, b2 in enumerate(_bases):
            for k, b3 in enumerate(_bases):
                CODON[b1 + b2 + b3] = _aas[(i * 16) + (j * 4) + k]

    def translate(dna):
        aa = []
        for i in range(0, len(dna) - 2, 3):
            aa.append(CODON.get(dna[i:i+3], "X"))
        return "".join(aa)

    n_win = 0
    n_hit = 0
    fam_hits = Counter()
    for chrom, lst in cds.items():
        seq = chrom_seqs.get(chrom, "")
        if not seq:
            continue
        for s, e, st in lst:
            cseq = seq[s-1:e] if st == "+" else seq[s-1:e].translate(comp)[::-1]
            prot = translate(cseq)
            for i in range(0, max(0, len(prot) - 59), 30):
                win = prot[i:i+60]
                if len(win) == 60 and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", win):
                    n_win += 1
                    hit = motif_hits(win)
                    if hit:
                        n_hit += 1
                        for f in hit:
                            fam_hits[f] += 1
    print(f"负集 (Yu11 CDS 蛋白 60aa 片段): {n_win} 窗口, motif 命中 {n_hit} "
          f"({n_hit/max(n_win,1)*100:.4f}%)")
    print(f"  负集逐家族命中: {dict(fam_hits)}")
    out.append({"set": "Yu11_CDS_60aa_neg", "n": n_win,
                "motif_hit_rate": n_hit})

    outp = ROOT / "results/08_benchmark/external_benchmark.tsv"
    with open(outp, "w") as f:
        f.write("set\tn\tlen_10_100\tcovered_by_any_sorf\t"
                "covered_by_ssp_candidate\tmotif_hit_rate\tnote\n")
        for r in out:
            f.write(f"{r['set']}\t{r['n']}\t{r.get('len_10_100','')}\t"
                    f"{r.get('covered_by_any_sorf','')}\t"
                    f"{r.get('covered_by_ssp_candidate','')}\t"
                    f"{r['motif_hit_rate']}\t{r.get('note','')}\n")
    print(f"\n✅ {outp}")


if __name__ == "__main__":
    main()
