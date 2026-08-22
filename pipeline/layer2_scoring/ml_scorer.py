"""
ml_scorer.py — LightGBM-based ML scoring for sORF coding potential

集成到PeptSesame Layer2评分系统的可替代证据通道。

架构:
  1. train.py → 训练LightGBM模型 → 保存为 .txt 文件
  2. predict.py → 加载模型 → 对sORF序列评分 → 0-1分
  3. 集成到 scoring_core.py → 作为Evidence F (ML证据) 加入加权评分

依赖: lightgbm, scikit-learn, numpy
"""
from __future__ import annotations

import json
import logging
import os
import pickle
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .feature_engineering import extract_all_features, extract_features_array, get_default_feature_names

logger = logging.getLogger(__name__)

# 默认特征名列表 (缓存)
_DEFAULT_FEATURE_NAMES: Optional[List[str]] = None


def _get_feature_names() -> List[str]:
    global _DEFAULT_FEATURE_NAMES
    if _DEFAULT_FEATURE_NAMES is None:
        _DEFAULT_FEATURE_NAMES = get_default_feature_names()
    return _DEFAULT_FEATURE_NAMES


# ===========================================================================
# LightGBM 预测器
# ===========================================================================

class LightGBMPredictor:
    """LightGBM模型预测器 — 对sORF序列给出0-1编码潜能得分。

    用法:
        predictor = LightGBMPredictor("models/lightgbm_sorf.txt")
        score = predictor.predict("MASKLCYFFLFLFLVLLSLPSSHCDD")
        # 返回 {"ml_score": 0.87, "probability": 0.92, "confidence": "high"}
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        model: Optional[object] = None,
        threshold_high: float = 0.7,
        threshold_medium: float = 0.4,
        feature_names: Optional[List[str]] = None,
    ):
        """
        Args:
            model_path: LightGBM模型文件路径 (.txt 格式)
            model: 可以直接传入已加载的模型对象 (优先于 model_path)
            threshold_high: HIGH置信度阈值 (default: 0.7)
            threshold_medium: MEDIUM置信度阈值 (default: 0.4)
            feature_names: 特征名列表 (默认自动获取)
        """
        self.threshold_high = threshold_high
        self.threshold_medium = threshold_medium
        self.feature_names = feature_names or _get_feature_names()

        if model is not None:
            self.model = model
            self.model_loaded = True
        elif model_path is not None and os.path.exists(model_path):
            self.model = self._load_model(model_path)
            self.model_loaded = True
        else:
            self.model = None
            self.model_loaded = False
            logger.warning("未找到LightGBM模型, ML评分降级为规则评分")

    def _load_model(self, model_path: str):
        """加载LightGBM模型。"""
        import lightgbm as lgb
        model = lgb.Booster(model_file=model_path)
        logger.info(f"LightGBM模型已加载: {model_path}")
        return model

    def is_ready(self) -> bool:
        """模型是否可用。"""
        return self.model_loaded and self.model is not None

    def predict(self, sequence: str) -> Dict:
        """对单条sORF序列预测编码潜能。

        Returns:
            Dict with keys:
                ml_score: float [0-1]
                probability: float [0-1] (模型原始概率)
                confidence: str ("high"|"medium"|"low")
        """
        if not self.is_ready():
            # 降级: 用规则评分
            return self._fallback_score(sequence)

        try:
            feats = extract_features_array(sequence, self.feature_names).reshape(1, -1)
            prob = self.model.predict(feats)[0]
            score = float(prob)  # LightGBM输出已经是0-1概率
        except Exception as e:
            logger.debug(f"ML评分失败, 降级到规则评分: {e}")
            return self._fallback_score(sequence)

        # 置信度分档
        if score >= self.threshold_high:
            confidence = "high"
        elif score >= self.threshold_medium:
            confidence = "medium"
        else:
            confidence = "low"

        return {
            "ml_score": round(score, 4),
            "probability": round(score, 4),
            "confidence": confidence,
        }

    def predict_batch(self, sequences: List[str]) -> List[Dict]:
        """批量预测。"""
        if not self.is_ready() or not sequences:
            return [self._fallback_score(s) for s in sequences]

        try:
            feats = np.array([
                extract_features_array(s, self.feature_names)
                for s in sequences
            ])
            probs = self.model.predict(feats)
            results = []
            for prob in probs:
                score = float(prob)
                if score >= self.threshold_high:
                    conf = "high"
                elif score >= self.threshold_medium:
                    conf = "medium"
                else:
                    conf = "low"
                results.append({
                    "ml_score": round(score, 4),
                    "probability": round(score, 4),
                    "confidence": conf,
                })
            return results
        except Exception as e:
            logger.warning(f"批量ML评分失败: {e}")
            return [self._fallback_score(s) for s in sequences]

    def _fallback_score(self, sequence: str) -> Dict:
        """降级方案: 基于简单的序列特征规则评分。"""
        from .scoring_core import SequenceFeatureComputer
        sfc = SequenceFeatureComputer()
        feat = sfc.compute_score(sequence)
        score = feat["combined"]
        seq_len = max(len(sequence), 1)

        # 调整: 太短或太长的降分
        if len(sequence) < 15:
            score *= 0.7
        elif len(sequence) > 80:
            score *= 0.8

        return {
            "ml_score": round(score, 4),
            "probability": round(score, 4),
            "confidence": "low",
        }


# ===========================================================================
# 训练数据准备
# ===========================================================================

def prepare_training_data(
    positive_fasta: str,
    negative_fasta: str,
    output_file: str,
    feature_names: Optional[List[str]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """从正负例FASTA准备训练数据。

    Args:
        positive_fasta: 正例sORF序列 (已知功能小肽)
        negative_fasta: 负例序列 (随机ORF/无表达区域)
        output_file: 特征矩阵输出路径 (.npz格式)
        feature_names: 特征名列表

    Returns:
        (X, y) 特征矩阵和标签
    """
    from Bio import SeqIO

    X_list, y_list = [], []
    fnames = feature_names or _get_feature_names()

    # 正例
    for r in SeqIO.parse(positive_fasta, "fasta"):
        seq = str(r.seq)
        feats = extract_features_array(seq, fnames)
        X_list.append(feats)
        y_list.append(1)

    # 负例
    for r in SeqIO.parse(negative_fasta, "fasta"):
        seq = str(r.seq)
        feats = extract_features_array(seq, fnames)
        X_list.append(feats)
        y_list.append(0)

    X = np.array(X_list)
    y = np.array(y_list)

    np.savez(output_file, X=X, y=y, feature_names=fnames)
    logger.info(f"训练数据已保存: {output_file} ({len(X)} samples)")
    return X, y


# ===========================================================================
# 训练入口
# ===========================================================================

def train_lightgbm(
    X: np.ndarray,
    y: np.ndarray,
    output_model: str = "models/lightgbm_sorf.txt",
    params: Optional[Dict] = None,
    feature_names: Optional[List[str]] = None,
) -> object:
    """训练LightGBM模型。

    Args:
        X: 特征矩阵 (n_samples, n_features)
        y: 标签 (0/1)
        output_model: 输出模型路径
        params: LightGBM参数
        feature_names: 特征名

    Returns:
        训练好的模型对象
    """
    import lightgbm as lgb

    default_params = {
        "objective": "binary",
        "metric": "auc",
        "boosting": "gbdt",
        "num_leaves": 31,
        "learning_rate": 0.05,
        "feature_fraction": 0.8,
        "bagging_fraction": 0.8,
        "bagging_freq": 5,
        "verbose": -1,
        "num_threads": 4,
        "min_data_in_leaf": 5,
        "is_unbalance": True,  # 正例通常远少于负例
    }
    if params:
        default_params.update(params)

    # 创建dataset
    fnames = feature_names or _get_feature_names()
    lgb_train = lgb.Dataset(X, y, feature_name=fnames)

    # 训练
    model = lgb.train(
        default_params,
        lgb_train,
        num_boost_round=500,
        callbacks=[lgb.log_evaluation(100)],
    )

    # 保存
    os.makedirs(os.path.dirname(output_model) or ".", exist_ok=True)
    model.save_model(output_model)

    # 保存参数
    param_path = output_model.replace(".txt", "_params.json")
    with open(param_path, "w") as f:
        json.dump(default_params, f, indent=2)

    logger.info(f"模型已保存: {output_model}")
    return model


def train_with_cv(
    X: np.ndarray,
    y: np.ndarray,
    output_dir: str = "models",
    n_folds: int = 5,
) -> Dict:
    """带交叉验证的训练 -- 用于评估。

    Returns:
        评估指标: {auc_mean, auc_std, precision, recall, f1}
    """
    import lightgbm as lgb
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score, precision_recall_fscore_support

    params = {
        "objective": "binary",
        "metric": "auc",
        "boosting": "gbdt",
        "num_leaves": 31,
        "learning_rate": 0.05,
        "verbose": -1,
        "is_unbalance": True,
    }
    fnames = _get_feature_names()

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    aucs, precisions, recalls, f1s = [], [], [], []

    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y)):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        lgb_train = lgb.Dataset(X_train, y_train, feature_name=fnames)
        lgb_val = lgb.Dataset(X_val, y_val, reference=lgb_train)

        model = lgb.train(
            params,
            lgb_train,
            num_boost_round=500,
            valid_sets=[lgb_val],
            callbacks=[lgb.early_stopping(20), lgb.log_evaluation(0)],
        )

        y_pred = (model.predict(X_val) > 0.5).astype(int)
        y_prob = model.predict(X_val)

        aucs.append(roc_auc_score(y_val, y_prob))
        p, r, f, _ = precision_recall_fscore_support(y_val, y_pred, average="binary")
        precisions.append(p)
        recalls.append(r)
        f1s.append(f)

        # 保存每个fold的模型
        fold_path = os.path.join(output_dir, f"lightgbm_fold{fold}.txt")
        model.save_model(fold_path)

    results = {
        "auc_mean": float(np.mean(aucs)),
        "auc_std": float(np.std(aucs)),
        "precision_mean": float(np.mean(precisions)),
        "recall_mean": float(np.mean(recalls)),
        "f1_mean": float(np.mean(f1s)),
        "n_folds": n_folds,
    }

    result_path = os.path.join(output_dir, "cv_results.json")
    with open(result_path, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"交叉验证结果: {results}")
    return results
