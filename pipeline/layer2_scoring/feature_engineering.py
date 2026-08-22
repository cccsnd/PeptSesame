"""
feature_engineering.py — sORF特征提取模块

将sORF肽段序列转化为数值特征向量，供LightGBM等ML模型使用。

特征分组:
  1. 氨基酸组成 (20维)
  2. 理化性质 (pI, MW, 疏水性, 电荷)
  3. 二肽频率 (400维, 可选PCA)
  4. 结构预测 (信号肽, TM域, 二硫键)
  5. 序列复杂度 (重复度, 熵)
"""
from __future__ import annotations

import math
import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np

# 氨基酸单字母码
AA_LETTERS = "ACDEFGHIKLMNPQRSTVWY"

# 氨基酸理化性质分组
AA_HYDROPHOBIC = {"A", "V", "I", "L", "M", "F", "W", "Y", "P"}
AA_POLAR = {"S", "T", "C", "N", "Q"}
AA_POSITIVE = {"K", "R", "H"}
AA_NEGATIVE = {"D", "E"}

# 分子量 (Da)
AA_MW = {
    "A": 89.09, "C": 121.15, "D": 133.10, "E": 147.13,
    "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.18,
    "K": 146.19, "L": 131.18, "M": 149.21, "N": 132.12,
    "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09,
    "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
}

# pKa值
AA_PKA_NH2 = 9.69  # 氨基端平均pKa
AA_PKA_COOH = 2.34  # 羧基端平均pKa
AA_PKA_SIDECHAIN = {
    "D": 3.86, "E": 4.25, "H": 6.00, "C": 8.33,
    "Y": 10.02, "K": 10.54, "R": 12.48,
}

# 疏水性标度 (Kyte-Doolittle)
AA_HYDRO_SCALE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "W": -0.9,
    "S": -0.8, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}


def _validate_aa(seq: str) -> str:
    """过滤非标准氨基酸，返回大写序列。"""
    return "".join(c for c in seq.upper() if c in AA_LETTERS)


# ==========================================================================
# 特征提取函数
# ==========================================================================

def aa_composition(seq: str, normalize: bool = True) -> Dict[str, float]:
    """氨基酸组成特征 (20维)。"""
    seq = _validate_aa(seq)
    if not seq:
        return {aa: 0.0 for aa in AA_LETTERS}
    counts = Counter(seq)
    total = len(seq) if normalize else 1.0
    return {aa: counts.get(aa, 0) / total for aa in AA_LETTERS}


def physicohemical_features(seq: str) -> Dict[str, float]:
    """理化性质特征 (10维)。"""
    seq = _validate_aa(seq)
    if not seq:
        return {"pI": 7.0, "MW": 0.0, "hydrophobicity": 0.0,
                "charge_at_pH7": 0.0, "aromacity": 0.0, "aliphatic_index": 0.0,
                "instability_index": 0.0, "cys_count": 0.0, "basic_ratio": 0.0,
                "acidic_ratio": 0.0}

    n = len(seq)
    # 分子量
    mw = sum(AA_MW.get(aa, 0) for aa in seq) - (n - 1) * 18.02  # 脱水
    # 平均疏水性
    hydro = sum(AA_HYDRO_SCALE.get(aa, 0) for aa in seq) / n
    # 净电荷 (pH 7.0)
    charge = 0.0
    for aa in seq:
        pka = AA_PKA_SIDECHAIN.get(aa, 7.0)
        if aa in AA_POSITIVE:
            charge += 1.0 / (1 + 10 ** (7.0 - pka))
        elif aa in AA_NEGATIVE:
            charge += -1.0 / (1 + 10 ** (pka - 7.0))
    # Cys数量
    cys_count = seq.count("C") / max(n, 1)
    # 碱性/酸性氨基酸比例
    basic = sum(1 for aa in seq if aa in "KRH") / max(n, 1)
    acidic = sum(1 for aa in seq if aa in "DE") / max(n, 1)
    # 芳香性 (Phe+Tyr+Trp)
    aroma = sum(1 for aa in seq if aa in "FYW") / max(n, 1)
    # 脂肪族指数
    aliphatic = sum(1 for aa in seq if aa in "AILV") / max(n, 1) * 100

    return {
        "pI": _estimate_pi(seq),
        "MW": round(mw, 2),
        "hydrophobicity": round(hydro, 4),
        "charge_at_pH7": round(charge, 4),
        "cys_ratio": round(cys_count, 4),
        "basic_ratio": round(basic, 4),
        "acidic_ratio": round(acidic, 4),
        "aromacity": round(aroma, 4),
        "aliphatic_index": round(aliphatic, 2),
        "length": n,
    }


def _estimate_pi(seq: str) -> float:
    """粗略估计等电点 (使用二分法)。"""
    def _net_charge(seq: str, pH: float) -> float:
        charge = 0.0
        # C端
        charge += -1.0 / (1 + 10 ** (pH - AA_PKA_COOH))
        # N端
        charge += 1.0 / (1 + 10 ** (AA_PKA_NH2 - pH))
        for aa in seq:
            pka = AA_PKA_SIDECHAIN.get(aa, 7.0)
            if aa in AA_POSITIVE:
                charge += 1.0 / (1 + 10 ** (pH - pka))
            elif aa in AA_NEGATIVE:
                charge += -1.0 / (1 + 10 ** (pka - pH))
        return charge

    lo, hi = 0.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if _net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def dipeptide_composition(seq: str) -> Dict[str, float]:
    """二肽频率特征 (400维)。"""
    seq = _validate_aa(seq)
    if not seq:
        return {a + b: 0.0 for a in AA_LETTERS for b in AA_LETTERS}
    dipeps = [seq[i:i+2] for i in range(len(seq) - 1)]
    counts = Counter(dipeps)
    total = len(dipeps)
    return {a + b: counts.get(a + b, 0) / max(total, 1)
            for a in AA_LETTERS for b in AA_LETTERS}


def sequence_complexity(seq: str) -> Dict[str, float]:
    """序列复杂度特征 (4维)。"""
    seq = _validate_aa(seq)
    if not seq:
        return {"shannon_entropy": 0.0, "repeat_ratio": 0.0,
                "low_complexity_score": 0.0, "n_glycosylation_sites": 0.0}

    n = len(seq)
    # Shannon熵
    counts = Counter(seq)
    entropy = -sum((c / n) * math.log2(c / n) for c in counts.values())
    # 重复比例 (相同氨基酸连续出现)
    repeats = sum(1 for i in range(1, n) if seq[i] == seq[i-1])
    # N-糖基化位点 (N-X-S/T)
    n_glyc = len(re.findall(r'N[^P][ST]', seq))
    return {
        "shannon_entropy": round(entropy / math.log2(20), 4) if n > 0 else 0.0,
        "repeat_ratio": round(repeats / max(n, 1), 4),
        "n_glycosylation_sites": n_glyc,
    }


def structural_prediction_features(seq: str) -> Dict[str, float]:
    """结构预测特征 (5维)。"""
    seq = _validate_aa(seq)
    if not seq:
        return {"signal_peptide_score": 0.0, "tm_helices": 0.0,
                "cysteine_rich": 0.0, "coiled_coil_prob": 0.0,
                "disorder_score": 0.0}

    n = len(seq)
    # 信号肽: N端疏水核心
    sp_score = 0.0
    if n >= 20:
        first20 = seq[:20]
        hydro_run = max((sum(1 for aa in first20[i:i+7] if aa in AA_HYDROPHOBIC)
                         for i in range(len(first20) - 6)), default=0)
        sp_score = min(hydro_run / 7.0, 1.0)

    # TM螺旋: 21aa窗口疏水
    tm_count = 0
    if n >= 21:
        i = 0
        while i <= n - 21:
            window = seq[i:i+21]
            if sum(1 for aa in window if aa in AA_HYDROPHOBIC) >= 14:
                tm_count += 1
                i += 21
            else:
                i += 1

    # Cys richness
    cys = seq.count("C")
    cys_rich = min(cys / 6.0 if cys >= 4 else 0.0, 1.0)

    return {
        "signal_peptide_score": round(sp_score, 4),
        "tm_helices": tm_count,
        "cysteine_rich_score": round(cys_rich, 4),
    }


# ==========================================================================
# 主特征提取函数
# ==========================================================================

def extract_all_features(seq: str, for_ml: bool = True) -> Dict[str, float]:
    """提取所有特征，返回展平的字典。

    Args:
        seq: 氨基酸序列
        for_ml: True=给ML模型用的完整特征(含二肽400维)
                False=给规则评分用的精简特征

    Returns:
        特征名字典
    """
    features = {}

    # 基础理化 (10维)
    features.update(physicohemical_features(seq))

    # 氨基酸组成 (20维)
    for k, v in aa_composition(seq).items():
        features[f"aa_{k}"] = v

    # 序列复杂度 (4维)
    features.update(sequence_complexity(seq))

    # 结构预测 (4维)
    features.update(structural_prediction_features(seq))

    # 二肽 (400维 — 仅ML模式)
    if for_ml:
        for k, v in dipeptide_composition(seq).items():
            features[f"dp_{k}"] = v

    return features


def extract_features_array(seq: str, feature_names: Optional[List[str]] = None) -> np.ndarray:
    """提取特征为numpy数组 (给sklearn/lightgbm直接使用)。"""
    features = extract_all_features(seq, for_ml=True)
    if feature_names is None:
        feature_names = sorted(features.keys())
    return np.array([features.get(f, 0.0) for f in feature_names], dtype=np.float32)


# ==========================================================================
# 特征名管理
# ==========================================================================

def get_default_feature_names() -> List[str]:
    """获取默认特征名列表 (用于保持一致的特征顺序)。"""
    return sorted(extract_all_features("MASKLCYFFLFLFLVLLSLPSSHCDD", for_ml=True).keys())


def feature_dimension() -> int:
    """返回总特征维度。"""
    return len(get_default_feature_names())
