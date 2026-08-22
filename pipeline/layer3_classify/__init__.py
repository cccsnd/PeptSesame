"""
PeptSesame — Layer3: Small Peptide Classification
==================================================

Classifies candidate sORF-encoded peptides into functional categories:
SSP families (CLE, RALF, CEP, etc.), miPEP, uORF, lncORF, and AMP.

Available components:

    SmallPeptideClassifier
        Main classification entry point — integrates all five sub-classifiers.

    SSPFamilyClassifier
        Motif-based matching against known SSP families.

    LocationClassifier
        GFF-based genomic context determination.

    AMPClassifier
        Antimicrobial peptide prediction (properties + database).

    ClassificationResult
        Dataclass holding classification results for one sORF.

    motif_profiles
        Curated regex patterns for SSP family detection
        (importable directly from pipeline.layer3_classify.motif_profiles).
"""

from pipeline.layer3_classify.classify import (
    SmallPeptideClassifier,
    SSPFamilyClassifier,
    LocationClassifier,
    AMPClassifier,
    ClassificationResult,
    ALL_CATEGORIES,
    CATEGORY_SSP,
    CATEGORY_MIPEP,
    CATEGORY_UORF,
    CATEGORY_LNCORF,
    CATEGORY_AMP,
)

__all__ = [
    "SmallPeptideClassifier",
    "SSPFamilyClassifier",
    "LocationClassifier",
    "AMPClassifier",
    "ClassificationResult",
    "ALL_CATEGORIES",
    "CATEGORY_SSP",
    "CATEGORY_MIPEP",
    "CATEGORY_UORF",
    "CATEGORY_LNCORF",
    "CATEGORY_AMP",
]
