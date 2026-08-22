#!/usr/bin/env python3
"""Layer3 classification (strict pipeline, motif + is_ida_like)

输入:
  Layer2: results_v2/02_layer2_scoring/<sp>/scored_sorfs.tsv (分数)
  Layer1: results_v2/01_layer1_sixframe/<sp>/sorfs.fa (序列)
输出:
  results_v2/03_layer3_classify/<sp>/classified_sorfs.tsv
自检: 行数 = Layer2 行数; is_ssp 计数; 家族计数; is_ida_like 计数

用法: python run_layer3.py <species>
"""
import os, sys, csv, re, time
from collections import Counter

ROOT = os.environ.get("PEPTSESAME_ROOT", ".")
R2 = f"{ROOT}/results_v2"

sys.path.insert(0, ROOT)
from pipeline.layer3_classify.motif_profiles import SSP_MOTIFS

FAMS = [f for f in SSP_MOTIFS if f != "IDA_LIKE"]  # 8 家族
IDA_LIKE_PAT = SSP_MOTIFS.get("IDA_LIKE")

def load_seqs(fa_path):
    """加载 FASTA → {id: seq} (id 无 sORF_ 前缀)"""
    seqs = {}
    cur = None
    buf = []
    with open(fa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur:
            seqs[cur] = "".join(buf)
    return seqs

def classify(species):
    src = f"{R2}/02_layer2_scoring/{species}/scored_sorfs.tsv"
    fa = f"{R2}/01_layer1_sixframe/{species}/sorfs.fa"
    if not os.path.exists(src) or not os.path.exists(fa):
        print(f"❌ {species} 输入缺失: {src} / {fa}")
        sys.exit(1)

    n_src = sum(1 for _ in open(src)) - 1
    print(f"=== {species} Layer3 classification (motif library) ===")
    print(f"Layer2 行数: {n_src:,}")

    print("加载序列...")
    seqs = load_seqs(fa)
    print(f"序列: {len(seqs):,}")

    out_dir = f"{R2}/03_layer3_classify/{species}"
    os.makedirs(out_dir, exist_ok=True)
    out = f"{out_dir}/classified_sorfs.tsv"

    t0 = time.time()
    fam_counts = Counter()
    ida_like_n = 0
    n_ssp = 0
    n_seq_missing = 0
    rows_out = 0
    with open(src) as fin, open(out, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        for f in FAMS:
            fieldnames.append(f"is_{f.lower()}")
        fieldnames += ["is_ida_like", "ssp_families", "is_ssp", "multi_class"]
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in reader:
            sid = r["seq_id"]
            seq = seqs.get(sid, "") or seqs.get(sid.replace("sORF_", ""), "")
            if not seq:
                n_seq_missing += 1
                # 保留行但不分类
                for f in FAMS:
                    r[f"is_{f.lower()}"] = "False"
                r["is_ida_like"] = "False"
                r["ssp_families"] = ""
                r["is_ssp"] = "False"
                r["multi_class"] = "False"
            else:
                fams = set()
                for f in FAMS:
                    if re.search(SSP_MOTIFS[f], seq):
                        fams.add(f)
                is_il = bool(re.search(IDA_LIKE_PAT, seq)) if IDA_LIKE_PAT else False
                for f in FAMS:
                    r[f"is_{f.lower()}"] = "True" if f in fams else "False"
                r["is_ida_like"] = "True" if is_il else "False"
                r["ssp_families"] = ";".join(sorted(fams)) if fams else ""
                r["is_ssp"] = "True" if fams else "False"
                r["multi_class"] = "True" if len(fams) > 1 else "False"
                fam_counts.update(fams)
                if is_il:
                    ida_like_n += 1
                if fams:
                    n_ssp += 1
            writer.writerow(r)
            rows_out += 1
            if rows_out % 200000 == 0:
                print(f"  ...{rows_out:,} ({time.time()-t0:.0f}s)", flush=True)
    print(f"\n完成: {rows_out:,} 行 ({time.time()-t0:.0f}s)")
    print(f"序列缺失: {n_seq_missing}")
    print(f"SSP 候选 (is_ssp=True): {n_ssp:,}")
    print(f"IDA_LIKE: {ida_like_n:,}")
    print(f"家族计数: {dict(fam_counts)}")

    # 自检
    n_out = sum(1 for _ in open(out)) - 1
    assert n_out == n_src, f"行数不一致: 输出 {n_out} vs 输入 {n_src}"
    print(f"✅ {species} Layer3 完成: {n_out:,} 行, 自检通过")

if __name__ == "__main__":
    classify(sys.argv[1])
