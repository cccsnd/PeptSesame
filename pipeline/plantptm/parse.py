"""PlantPTM result parsing — batch JSON -> site rows + sequon sanity.

Consumption layer for raw batch results saved by PlantPTMClient workflows:
    load_raw_rows(rawdir)      -> flattened site rows from batch_*.json
    sanity_ngly_sequon(rows, seqs) -> per-confidence-tier N-X-S/T consistency
                                     (high-conf tiers ~100% consistent, Low is
                                     noise — filter >= High before interpretation)
"""
from __future__ import annotations

import json
from pathlib import Path

# 9类PTM → 中文/含义 (for reports) and modified residues (for sanity checks)
PTM_CN = {
    "Ngly": "N-糖基化", "Sacy": "S-酰化", "Khib": "2-羟基异丁酰化",
    "Kcr": "巴豆酰化", "Ksucc": "琥珀酰化", "Kmal": "丙二酰化",
    "Kac": "乙酰化", "Kub": "泛素化", "pho": "磷酸化",
}
PTM_RES = {"Ngly": "N", "Sacy": "C", "Khib": "K", "Kcr": "K", "Ksucc": "K",
           "Kmal": "K", "Kac": "K", "Kub": "K", "pho": "STY"}


def load_raw_rows(rawdir: Path) -> list[dict]:
    """Flatten all batch_*.json in rawdir into site-level rows."""
    rows = []
    for p in sorted(rawdir.glob("batch_*.json")):
        with open(p) as f:
            result = json.load(f)
        rows.extend(json.loads(result.get("full_data_json") or "[]"))
    return rows


def sanity_ngly_sequon(rows: list[dict], seqs: dict) -> dict:
    """Ngly sites vs the N-X-S/T (X!=P) sequon, stratified by confidence tier.

    Model-vs-biochemistry cross-validation: an aggregate rate hides that
    high-confidence calls are reliable while Low is pure noise.
    """
    ngly = [r for r in rows if r["PTM type"] == "Ngly"]
    levels = ["Extremely high", "High", "Medium", "Low", "Non-PTM"]
    per_level = {lv: {"total": 0, "ok": 0, "bad": 0, "no_ctx": 0, "rate": None}
                 for lv in levels}
    examples_bad = []
    for r in ngly:
        seq = seqs.get(r["Protein"])
        lv = r.get("Confidence", "Non-PTM")
        per_level[lv]["total"] += 1
        if not seq:
            per_level[lv]["no_ctx"] += 1
            continue
        pos = int(r["Position"]) - 1  # 0-based
        if pos + 2 >= len(seq):
            per_level[lv]["no_ctx"] += 1
            continue
        x, st = seq[pos + 1], seq[pos + 2]
        if x == "P":
            per_level[lv]["bad"] += 1
            if len(examples_bad) < 5:
                examples_bad.append((r["Protein"], r["Position"], r["Confidence"],
                                     seq[max(0, pos - 3):pos + 7]))
        elif st in "ST":
            per_level[lv]["ok"] += 1
        else:
            per_level[lv]["bad"] += 1
            if len(examples_bad) < 5:
                examples_bad.append((r["Protein"], r["Position"], r["Confidence"],
                                     seq[max(0, pos - 3):pos + 7]))
    for lv in levels:
        d = per_level[lv]
        if (d["ok"] + d["bad"]) > 0:
            d["rate"] = d["ok"] / (d["ok"] + d["bad"]) * 100
        else:
            d["rate"] = None
    return {"total": len(ngly), "per_level": per_level, "examples_bad": examples_bad}
