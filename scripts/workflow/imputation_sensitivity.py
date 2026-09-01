import os
#!/usr/bin/env python3
"""0.5 neutral imputation sensitivity (R1-4/评估): 三种缺失值处理方案对比

A. 当前方案: 缺失通道 = 0.5 (neutral)
B. 可用通道重归一化: 缺失通道不计入分子分母 (C = Σ_avail w_i S_i / Σ_avail w_i)
C. 缺失指示: 缺失通道直接剔除 (权重重分配, 等价于 B 的确定性版本)

对 Yu11 SSP 候选: 比较 expression 通道有真实值 (当前) vs 模拟缺失 (0.5 / 剔除)
时排序的稳定性 (Spearman ρ / top-k Jaccard)。
"""
import csv, math
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
V5 = ROOT / "results"
OUT = V5 / "08_benchmark/imputation_sensitivity_v5.tsv"

# 5 通道归一化权重 (去 ML 后)
W = {"sequence_features": 0.20/0.75, "conservation": 0.15/0.75,
     "expression": 0.10/0.75, "structural": 0.15/0.75, "motif": 0.15/0.75}
CH = ["sequence_features", "conservation", "expression", "structural", "motif"]


def load_ssp():
    rows = []
    with open(V5 / "03_layer3_classify/Yu11/classified_sorfs.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["is_ssp"] == "True":
                rows.append({c: float(r[c]) for c in CH})
    return rows


def score_full(r):
    return sum(W[c] * r[c] for c in CH)


def score_no_expr(r):  # expression 缺失 → 0.5
    rr = dict(r); rr["expression"] = 0.5
    return sum(W[c] * rr[c] for c in CH)


def score_renorm(r):  # expression 缺失 → 剔除该通道, 重归一化
    wsum = sum(W[c] for c in CH if c != "expression")
    return sum(W[c] * r[c] for c in CH if c != "expression") / wsum


def rho(a, b):
    n = len(a)
    ra = [0]*n
    for k, i in enumerate(sorted(range(n), key=lambda i: a[i], reverse=True)):
        ra[i] = k
    rb = [0]*n
    for k, i in enumerate(sorted(range(n), key=lambda i: b[i], reverse=True)):
        rb[i] = k
    ma, mb = sum(ra)/n, sum(rb)/n
    cov = sum((ra[i]-ma)*(rb[i]-mb) for i in range(n))
    va = math.sqrt(sum((x-ma)**2 for x in ra))
    vb = math.sqrt(sum((x-mb)**2 for x in rb))
    return cov/(va*vb) if va and vb else 0.0


def jac(a, b, k):
    ta = set(sorted(range(len(a)), key=lambda i: a[i], reverse=True)[:k])
    tb = set(sorted(range(len(b)), key=lambda i: b[i], reverse=True)[:k])
    return len(ta & tb)/len(ta | tb)


def main():
    rows = load_ssp()
    n = len(rows)
    print(f"Yu11 SSP 候选: {n}")
    sA = [score_full(r) for r in rows]      # 当前 (expression 有值)
    sB = [score_no_expr(r) for r in rows]   # expression 模拟缺失=0.5
    sC = [score_renorm(r) for r in rows]    # expression 剔除重归一化

    out = [("comparison", "spearman_rho", "top50_jaccard", "top100_jaccard")]
    for name, sa, sb in [("A_vs_B(0.5)", sA, sB), ("A_vs_C(renorm)", sA, sC), ("B_vs_C", sB, sC)]:
        r = rho(sa, sb)
        j50 = jac(sa, sb, 50)
        j100 = jac(sa, sb, 100)
        print(f"  {name}: ρ={r:.4f}  top50 J={j50:.3f}  top100 J={j100:.3f}")
        out.append((name, f"{r:.4f}", f"{j50:.3f}", f"{j100:.3f}"))

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerows(out)
    print(f"✅ {OUT}")


if __name__ == "__main__":
    main()
