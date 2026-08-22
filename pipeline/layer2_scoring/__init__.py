"""
PeptSesame — Layer2: Multi-Evidence Scoring System
===================================================

Provides the EvidenceScorer class that integrates six computational
evidence channels to assess sORF coding potential and biological
relevance, replacing experimental evidence (Ribo-seq / MS).

Evidence channels:
    A. Sequence features (k-mer composition, GC bias, hexamer scores)
    B. Cross-species conservation (BLASTP against related species peptides)
    C. RNA-seq expression evidence (if BAM/alignment files available)
    D. Structural features (signal peptide, transmembrane, cysteine-richness)
    E. Known motif / domain detection
    F. ML-based coding potential (LightGBM — if model is available)

Available components:

    EvidenceScorer
        Main class combining all six evidence channels.

    EvidenceScores
        Dataclass holding per-channel and aggregated scores.

    SequenceFeatureComputer
        Evidence A — k-mer composition, GC bias, hexamer scores.

    ConservationComputer
        Evidence B — cross-species BLASTP conservation.

    ExpressionComputer
        Evidence C — RNA-seq expression (RPKM/TPM).

    StructuralFeatureComputer
        Evidence D — signal peptide, TM domains, CRP.

    MotifComputer
        Evidence E — known motif / domain matching.

    MlScorerComputer
        Evidence F — ML-based coding potential (LightGBM).

    LightGBMPredictor
        ML prediction engine (from ml_scorer module).
"""

from pipeline.layer2_scoring.scoring_core import (
    EvidenceScorer,
    EvidenceScores,
    SequenceFeatureComputer,
    ConservationComputer,
    ExpressionComputer,
    StructuralFeatureComputer,
    MotifComputer,
    MlScorerComputer,
    DEFAULT_WEIGHTS,
)

from pipeline.layer2_scoring.ml_scorer import (
    LightGBMPredictor,
    prepare_training_data,
    train_lightgbm,
    train_with_cv,
)

from pipeline.layer2_scoring.feature_engineering import (
    extract_all_features,
    extract_features_array,
    get_default_feature_names,
    feature_dimension,
    aa_composition,
    physicohemical_features,
    dipeptide_composition,
)

__all__ = [
    "EvidenceScorer",
    "EvidenceScores",
    "SequenceFeatureComputer",
    "ConservationComputer",
    "ExpressionComputer",
    "StructuralFeatureComputer",
    "MotifComputer",
    "MlScorerComputer",
    "LightGBMPredictor",
    "DEFAULT_WEIGHTS",
    "prepare_training_data",
    "train_lightgbm",
    "train_with_cv",
    "extract_all_features",
    "extract_features_array",
    "get_default_feature_names",
    "feature_dimension",
    "aa_composition",
    "physicohemical_features",
    "dipeptide_composition",
]
