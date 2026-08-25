#!/usr/bin/env python3
"""Motif-hit empirical null over the Yu11 search space (1,000 permutations).

Reproduces and extends the motif-hit background analysis to 1,000 permutations:
- 500,000 sORFs sampled without replacement from the Yu11 Layer 1 library (seed 42)
- per-sequence amino-acid shuffling (preserves composition and length)
- Layer 3 motif pipeline applied to each shuffled set
- per-family hit rates normalized to the full library: (hits/500,000) x 4,451,664
- empirical P = (1 + null >= observed) / (1 + 1,000), one-sided
- Benjamini-Hochberg correction across the nine motif frameworks

Observed full-library hit counts are read from the original analysis
(real_hits column of empirical_null_v20260818.tsv).

Run: python scripts/workflow/null_model_motif.py [--perms 1000] [--workers 32]
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import random
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pipeline.layer3_classify.motif_profiles import SSP_MOTIFS  # noqa: E402

COMPILED = {fam: re.compile(p) for fam, p in SSP_MOTIFS.items()}

SAMPLE_N = 500_000
FULL_LIB = 4_451_664
SEED = 42
ORIGINAL = ROOT / "results/08_benchmark/empirical_null_v20260818.tsv"
OUT = ROOT / "results/08_benchmark/empirical_null_motif_1000x.tsv"
SORFS_FA = ROOT / "results/01_layer1_sixframe/Yu11/sorfs.fa"


def load_sequences() -> list[str]:
    seqs: list[str] = []
    buf: list[str] = []
    for line in SORFS_FA.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf).upper())
                buf = []
            continue
        if line.strip():
            buf.append(line.strip())
    if buf:
        seqs.append("".join(buf).upper())
    return seqs


def shuffled_hits(seqs: list[str], rng: random.Random) -> dict[str, int]:
    hits = {fam: 0 for fam in COMPILED}
    for seq in seqs:
        aa = list(seq)
        rng.shuffle(aa)
        shuffled = "".join(aa)
        for fam, pat in COMPILED.items():
            if pat.search(shuffled):
                hits[fam] += 1
    return hits


def worker(task: tuple[int, list[str], int]) -> dict[str, int]:
    perm_idx, sample, seed = task
    rng = random.Random(seed + perm_idx)
    return shuffled_hits(sample, rng)


def bh_correction(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, pvals[i] * n / (rank + 1))
        q[i] = running
    return q


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    seqs = load_sequences()
    print(f"库: {len(seqs):,} 序列")
    rng = random.Random(args.seed)
    sample = rng.sample(seqs, SAMPLE_N)
    del seqs
    print(f"固定样本: {len(sample):,} (seed {args.seed})")

    real = {}
    for line in ORIGINAL.read_text(encoding="utf-8").splitlines()[1:]:
        fam, hits, *_ = line.split("\t")
        real[fam] = int(hits)

    tasks = [(i, sample, args.seed) for i in range(args.perms)]
    with mp.Pool(args.workers) as pool:
        results = pool.map(worker, tasks, chunksize=8)

    fams = list(COMPILED.keys())
    null_counts = {fam: [] for fam in fams}
    for res in results:
        for fam in fams:
            null_counts[fam].append(res[fam])

    rows = []
    for fam in fams:
        obs = real.get(fam, 0)
        nulls = sorted(null_counts[fam])
        mean = sum(nulls) / len(nulls)
        lo = nulls[int(0.025 * len(nulls))]
        hi = nulls[int(0.975 * len(nulls)) - 1]
        # normalize sampled counts to full library
        scale = FULL_LIB / SAMPLE_N
        norm = [n * scale for n in nulls]
        nmean = sum(norm) / len(norm)
        nlo = norm[int(0.025 * len(norm))]
        nhi = norm[int(0.975 * len(norm)) - 1]
        p = (1 + sum(1 for n in norm if n >= obs)) / (1 + len(norm))
        fold = obs / nmean if nmean > 0 else float("inf")
        rows.append((fam, obs, nmean, nlo, nhi, p, 0.0, fold))

    qvals = bh_correction([r[5] for r in rows])
    rows = [(f, o, m, lo, hi, p, q, fold) for (f, o, m, lo, hi, p, _, fold), q in zip(rows, qvals)]

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("family\treal_hits\tnull_mean\tnull_ci_low\tnull_ci_high\tempirical_p\tbh_q\tfold_enrichment\n")
        for fam, obs, m, lo, hi, p, q, fold in rows:
            fh.write(f"{fam}\t{obs}\t{m:.1f}\t{lo:.1f}\t{hi:.1f}\t{p:.4f}\t{q:.4f}\t{fold:.2f}\n")
    print(f"✅ {OUT}")
    for r in rows:
        print(f"  {r[0]:9s} obs={r[1]:6d} null_mean={r[2]:7.1f} P={r[5]:.4f} q={r[6]:.4f} fold={r[7]:.2f}")


if __name__ == "__main__":
    main()
