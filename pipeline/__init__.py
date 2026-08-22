"""
PeptSesame Pipeline
===================
Plant small peptide (small ORF) mining pipeline.

Layers:
    Layer 1: Six-frame translation + ORF filtering (genome-level sORF discovery)
    Layer 2: Multi-evidence scoring system (coding potential, conservation, expression)
    Layer 3: Small peptide classification (SSP families, miPEP, uORF, lncORF, AMP)
    Layer 4: Functional prediction (GO/KEGG, stress response, receptor pairing)
"""

__version__ = "0.1.0"
