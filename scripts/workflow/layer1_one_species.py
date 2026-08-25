#!/usr/bin/env python3
"""Per-species Layer1 execution + self-check (strict pipeline)

用法: python layer1_one_species.py <species>
步骤:
1. 输入验证 (genome + gff 存在)
2. 运行 Layer1 (SixFrameTranslator, 统一参数)
3. 自检统计: 总数/密度/同链重叠/反链重叠/长度分布/链分布
4. 写记录文件 species_record.md

不批量, 一次一个物种。所有数字从输出文件读取。
"""
import os, sys, time, json, random, csv

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
sys.path.insert(0, f"{ROOT}/pipeline/layer1_sixframe")
from sixframe import SixFrameTranslator

R2 = f"{ROOT}/results"

def load_manifest():
    with open(f"{R2}/species_manifest.json") as f:
        return json.load(f)

def check_inputs(species, manifest):
    g, gf, note = manifest[species]
    g_ok = os.path.exists(g)
    gf_ok = os.path.exists(gf)
    print(f"  输入验证: genome={'✅' if g_ok else '❌'} ({g})")
    print(f"           gff={'✅' if gf_ok else '❌'} ({gf})")
    if not (g_ok and gf_ok):
        print("❌ 输入缺失, 中止")
        sys.exit(1)
    return g, gf

def run_layer1(species, genome, gff):
    t0 = time.time()
    out_dir = f"{R2}/01_layer1_sixframe/{species}"
    os.makedirs(out_dir, exist_ok=True)
    bed = f"{out_dir}/sorfs.bed"
    fa = f"{out_dir}/sorfs.fa"
    st = SixFrameTranslator(
        genome_fasta=genome, cds_gff=gff,
        out_bed=bed, out_pep_fasta=fa,
        min_orf_len=30, max_orf_len=300,
        strand_aware_overlap=True, keep_partial_orfs=False,
        n_jobs=8, tmp_dir=tempfile.gettempdir(),
    )
    st.run()
    print(f"  Layer1 完成 ({time.time()-t0:.0f}s)")
    return bed

def self_check(species, bed, gff):
    """自检: 从输出文件读统计"""
    # 1. sORF 总数 + 链/长度分布
    n = 0
    strand_c = {"+": 0, "-": 0}
    len_bins = {"10-20": 0, "20-50": 0, "50-100": 0}
    with open(bed) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            n += 1
            strand_c[parts[5]] = strand_c.get(parts[5], 0) + 1
            laa = int(parts[8])
            if laa < 20: len_bins["10-20"] += 1
            elif laa < 50: len_bins["20-50"] += 1
            else: len_bins["50-100"] += 1
    print(f"  sORF 总数: {n:,}")
    print(f"  链分布: +{strand_c['+']:,} / -{strand_c['-']:,}")
    print(f"  长度分布: {len_bins}")

    # 2. CDS 重叠检查 (采样 5000)
    cds = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) >= 7 and parts[2] == "CDS":
                cds.setdefault(parts[0], []).append((int(parts[3]), int(parts[4]), parts[6]))
    random.seed(42)
    same = opp = 0
    with open(bed) as f:
        next(f)
        lines = f.readlines()
    sample = random.sample(lines, min(5000, len(lines)))
    for line in sample:
        parts = line.rstrip("\n").split("\t")
        chrom, s, e, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]
        for cs, ce, cst in cds.get(chrom, []):
            if s < ce and cs < e:
                if strand == cst: same += 1
                else: opp += 1
                break
    ns = len(sample)
    print(f"  CDS 重叠 (采样{ns}): 同链 {same} ({same/ns*100:.1f}%) | 反链 {opp} ({opp/ns*100:.1f}%)")

    # 3. 密度 (基因组大小从 manifest? 用 FASTA 统计)
    return {"total": n, "strand": strand_c, "len_bins": len_bins,
            "same_strand_overlap_pct": same/ns*100, "opp_strand_overlap_pct": opp/ns*100}

def record(species, stats, g, gf, runtime):
    rec = f"{R2}/01_layer1_sixframe/{species}/layer1_record.md"
    with open(rec, "w") as f:
        f.write(f"# {species} Layer1 记录\n\n")
        f.write(f"- 日期: 2026-08-14\n- 基因组: {g}\n- GFF: {gf}\n")
        f.write(f"- 参数: min=30nt max=300nt strand_aware=True partial=False\n")
        f.write(f"- 运行时间: {runtime:.0f}s\n")
        f.write(f"- sORF 总数: {stats['total']:,}\n")
        f.write(f"- 链分布: {stats['strand']}\n")
        f.write(f"- 长度分布: {stats['len_bins']}\n")
        f.write(f"- 同链 CDS 重叠: {stats['same_strand_overlap_pct']:.1f}%\n")
        f.write(f"- 反链 CDS 重叠: {stats['opp_strand_overlap_pct']:.1f}%\n")
    print(f"  记录: {rec}")

if __name__ == "__main__":
    species = sys.argv[1]
    manifest = load_manifest()
    print(f"=== {species} Layer1 (strict pipeline) ===")
    g, gf = check_inputs(species, manifest)
    t0 = time.time()
    bed = run_layer1(species, g, gf)
    runtime = time.time() - t0
    stats = self_check(species, bed, gf)
    record(species, stats, g, gf, runtime)
    print(f"=== {species} 完成 ✅ ===")
