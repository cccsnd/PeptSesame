#!/usr/bin/env python3
"""null_model_100x.py — sORF 密度 null model, 每条染色体前 5Mb 子集 × 100 次置换

Method:
- 对 Yu11 / Arabidopsis / maize 每条染色体的前 5 Mb:
  * 保留单核苷酸组成洗牌 (二核苷酸保留更稳但成本高; mononucleotide-shuffle baseline)
  * 六框扫描 30-300nt ORF 计数 (简化计数法: 正向 3 框 + 反向 3 框, ATG 起始 + 同框终止)
  * 100 次置换 → 均值 + 95% CI (2.5%/97.5% 分位)
- 输出: 观测密度 vs null 密度 (mean, 95% CI) + 观测/null 比值 (mean, CI)

输出: results_v2/08_benchmark/null_model_100x.tsv
"""
import random
import os
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
R2 = ROOT / "results_v2"
OUT = R2 / "08_benchmark"

SPECIES = {
    "Yu11": "results_v2/00_inputs/genomes/Yu11.fasta",
    "Arabidopsis": "results_v2/00_inputs/genomes/Arabidopsis.fasta",
    "Maize": "results_v2/00_inputs/genomes/Maize_B73.fasta",
}
N_PERM = 100
SUB = 5_000_000  # 每条染色体前 5 Mb


def load_chroms(fasta: Path, max_bp: int):
    """返回 [(chrom, seq[:max_bp])]"""
    out = []
    cur, buf = None, []
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur and buf:
                    out.append((cur, "".join(buf)[:max_bp]))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if cur and buf:
        out.append((cur, "".join(buf)[:max_bp]))
    return out


def count_sorfs(seq, min_len=30, max_len=300):
    """六框 sORF 计数 (find 跳跃, 快): 每框找 ATG + 同框终止, 长度 min~max"""
    comp = str.maketrans("ACGT", "TGCA")
    rev = seq.translate(comp)[::-1]
    n = 0
    stops = ("TAA", "TAG", "TGA")
    for strand in (seq, rev):
        L = len(strand)
        for frame in range(3):
            i = frame
            while i < L - 2:
                atg = strand.find("ATG", i)
                if atg < 0:
                    break
                # 找同框终止
                j = atg + 3
                while j < L - 2:
                    codon = strand[j:j+3]
                    if codon in stops:
                        ln = j + 3 - atg
                        if min_len <= ln <= max_len:
                            n += 1
                        break
                    j += 3
                i = j if j < L - 2 else atg + 1
    return n


def shuffle_keep_nt(seq, rng):
    """保留单核苷酸组成的洗牌"""
    s = list(seq)
    rng.shuffle(s)
    return "".join(s)


def main():
    rng = random.Random(42)
    with open(OUT / "null_model_100x.tsv", "w") as f:
        f.write("species\tchrom\tobs_density_per_mb\tnull_mean_density_per_mb\t"
                "null_ci_low\tnull_ci_high\tratio_mean\tratio_ci_low\tratio_ci_high\n")
        f.flush()
        for sp, fa_path in SPECIES.items():
            fa = ROOT / fa_path
            if not fa.exists():
                print(f"⚠️ 缺失 {fa}, 跳过 {sp}")
                continue
            print(f"=== {sp}: 加载染色体 (前 {SUB/1e6:.0f} Mb) ===", flush=True)
            chroms = [s for _, s in load_chroms(fa, SUB) if len(s) >= 100_000]
            t0 = time.time()
            for ci, seq in enumerate(chroms, 1):
                obs = count_sorfs(seq)
                obs_dens = obs / (len(seq) / 1e6)
                nulls = []
                for k in range(N_PERM):
                    sh = shuffle_keep_nt(seq, rng)
                    nulls.append(count_sorfs(sh) / (len(seq) / 1e6))
                nulls.sort()
                nm = sum(nulls) / len(nulls)
                lo, hi = nulls[int(0.025 * N_PERM)], nulls[int(0.975 * N_PERM) - 1]
                ratios = [obs_dens / x for x in nulls]
                ratios.sort()
                rm = sum(ratios) / len(ratios)
                rlo, rhi = ratios[int(0.025 * N_PERM)], ratios[int(0.975 * N_PERM) - 1]
                f.write(f"{sp}\t{ci}\t{obs_dens:.1f}\t{nm:.1f}\t{lo:.1f}\t{hi:.1f}\t"
                        f"{rm:.3f}\t{rlo:.3f}\t{rhi:.3f}\n")
                f.flush()
                print(f"  chr{ci}: obs={obs_dens:.0f}/Mb null={nm:.0f}/Mb "
                      f"ratio={rm:.3f} [{rlo:.3f},{rhi:.3f}] "
                      f"({time.time()-t0:.0f}s)", flush=True)
    print("✅ null_model_100x.tsv 完成")


if __name__ == "__main__":
    main()
