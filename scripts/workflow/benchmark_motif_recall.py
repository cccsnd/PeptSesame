#!/usr/bin/env python3
"""benchmark_motif_recall.py — SSP 基序 recall/FPR/PPV 基准 (08_benchmark, current pipeline)

Method (identical to the earlier gen_tableS6_motif_benchmark.py):
- 金标准: 拟南芥已知 SSP 家族成员序列 (CLE 6 / RALF 1 / PSK 5 / IDA 1)
- 匹配: motif library (pipeline/layer3_classify/motif_profiles.py SSP_MOTIFS),
  与 run_layer3.py 生产逻辑一致 (re.search, 大小写敏感与生产一致)
- recall = 金标准命中比例; FPR = 5,000 条随机 60aa 肽命中比例
- PPV 估计 = recall / (recall + fp_rate×99) (假设金标准正例占候选池 1%)

输出: results/08_benchmark/TableS6_motif_benchmark.tsv
"""
import random
import re
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("PEPTSESAME_ROOT", "."))
sys.path.insert(0, str(ROOT))
from pipeline.layer3_classify.motif_profiles import SSP_MOTIFS

OUT_DIR = ROOT / "results/08_benchmark"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# gold standard (as in the earlier gen_tableS6_motif_benchmark.py; Arabidopsis known family members)
GOLD = {
    "CLE": ["MANLKFLLCLFLICVSLSRSSASRPMFPNADGIKRGRMMIEAEEVLKASMEKLMERGFNESMRLSPGGPDPRHH",
            "MAKLSFTFCFLLFLLLSSIAAGSRPLEGARVGVKVRGLSPSIEATSPTVEDDQAAGSHGKSPERLSPGGPDPQHH",
            "MASLKLWVCLVLLLVLELTSVHECRPLVAEERFSGSSRLKKIRRELFERLKEMKGRSEGEETILGNTLDSKRLSPGGPDPRHH",
            "MATLILKQTLIILLIIFSLQTLSSQARILRSYRAVSMGNMDSQVLLHELGFDLSKFKGHNERRFLVSSDRVSPGGPDPQHH",
            "MDSKSFLLLLLLFCFLFLHDASDLTQAHAHVQGLSNRKMMMMKMESEWVGANGEAEKAKTKGLGLHEELRTVPSGPDPLHH",
            "MAAMKYKGSVFIILVILLLSSSLLAHSSSTKSFFWLGETQDTKAMKKEKKIDGGTANEVEERQVPTGSDPLHHKHIPFTP"],
    "RALF": ["MASKLCYFFLFLFLVLLSLPSSHATCNLKDCVNEADASNLTAMRAVSVPVSVSKGLGDEELTQSVYVSCVDGASPKRVPCNRRGFANVPRYISY"],
    "PSK": ["MKTKSEVLIFFFTLVLLLSMASSVILREDGFAPPKPSPTTHEKASTKGDRDGVECKNSDSEEECLVKKTVAAHTDYIYTQDLNLSP",
            "MANVSALLTIALLLCSTLMCTARPEPAISISITTAADPCNMEKKIEGKLDDMHMVDENCGADDEDCLMRRTLVAHTDYIYTQKKKHP",
            "MKQSLCLAVLFLILSTSSSAIRRGKEDQEINPLVSATSVEEDSVNKLMGMEYCGEGDEECLRRRMMTESHLDYIYTQHHKH",
            "MGKFTTIFIMALLLCSTLTYAARLTPTTTTALSRENSVKEIEGDKVEEESCNGIGEEECLIRRSLVLHTDYIYTQNHKP",
            "MVKFTTFLCIIALLLCSTLTHASARLNPTSVYPEENSFKKLEQGEVICEGVGEEECFLIRRTLVAHTDYIYTQNHNP"],
    "IDA": ["MRNNHSLRLQLWFRTLFTVGVVTLMIDAFVLQNNKEDDKTKEITTAVNMNNSDAKEIQQELEDGSRNDDLSYVASKRKVPRGPDPIHNRRAGNSRRPPGRA"],
}

random.seed(42)
AA = "ACDEFGHIKLMNPQRSTVWY"
N_RAND = 5000

rows = []
print("family\tn_gold\trecall\ttp_rate_rand\tppv_est")
for fam, seqs in GOLD.items():
    pat = SSP_MOTIFS[fam]
    rec = sum(1 for s in seqs if re.search(pat, s))
    recall = rec / len(seqs)
    rand_hits = 0
    for _ in range(N_RAND):
        pep = "".join(random.choice(AA) for _ in range(60))
        if re.search(pat, pep):
            rand_hits += 1
    fp_rate = rand_hits / N_RAND
    ppv = recall / (recall + fp_rate * 99) if (recall + fp_rate * 99) > 0 else 0
    rows.append((fam, len(seqs), recall, fp_rate, ppv))
    print(f"{fam}\t{len(seqs)}\t{recall:.2f}\t{fp_rate:.4f}\t{ppv:.3f}")

out = OUT_DIR / "TableS6_motif_benchmark.tsv"
with open(out, "w") as f:
    f.write("family\tn_gold_standard\trecall\ttp_rate_random60aa\tppv_estimate\n")
    for fam, n, rec, fp, ppv in rows:
        f.write(f"{fam}\t{n}\t{rec:.2f}\t{fp:.4f}\t{ppv:.3f}\n")
print(f"\n✅ {out}")
