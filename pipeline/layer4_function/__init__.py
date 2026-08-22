"""
PeptSesame — Layer4: Functional Prediction & Visualization
===========================================================

Adds functional context to sORF-encoded peptides identified by Layers 1-3:

    1. GO annotation       — homology-based Gene Ontology term transfer
    2. KEGG pathway mapping — BLASTP-based KEGG pathway assignment
    3. Stress association   — RNA-seq differential expression integration
    4. Receptor pairing     — LRR-RLK partner prediction for SSPs (CLE, RALF, CEP)
    5. Tissue specificity   — multi-tissue expression pattern analysis
    6. Conservation profile — cross-species conservation heatmap
    7. Visual summary       — matplotlib summary plots (histograms, pies, heatmaps)

Typical usage::

    from pipeline.layer4_function.functional import FunctionalAnnotator

    annotator = FunctionalAnnotator()
    go_terms = annotator.go_annotation("MASSKLCYFFLFLFLV...")
    pathways = annotator.kegg_mapping("MASSKLCYFFLFLFLV...")
    plots = annotator.generate_plots(sorf_annotations, output_dir="plots/")

Available components:

    FunctionalAnnotator
        Main class for functional prediction and visualisation.

    GOAnnotation
        Dataclass for a single GO term assignment.

    KEGGResult
        Dataclass for a KEGG pathway assignment.

    ReceptorPairingResult
        Dataclass for a predicted receptor-peptide pair.

    StressAssociationResult
        Dataclass for stress-response association details.
"""

from pipeline.layer4_function.functional import (
    FunctionalAnnotator,
    GOAnnotation,
    KEGGResult,
    ReceptorPairingResult,
    StressAssociationResult,
    BASE_GO_TERMS,
    PLANT_GO_TERMS,
    KEGG_PLANT_PATHWAYS,
    SSP_RECEPTOR_PAIRS,
)

__all__ = [
    "FunctionalAnnotator",
    "GOAnnotation",
    "KEGGResult",
    "ReceptorPairingResult",
    "StressAssociationResult",
    "BASE_GO_TERMS",
    "PLANT_GO_TERMS",
    "KEGG_PLANT_PATHWAYS",
    "SSP_RECEPTOR_PAIRS",
]
