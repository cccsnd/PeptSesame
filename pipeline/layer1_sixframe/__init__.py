"""
PeptSesame Pipeline — Layer 1: Six-Frame Genome Translation & sORF Filtering.

This module performs six-frame translation of a genome FASTA, identifies
small open reading frames (sORFs) between 30-300 nt (10-100 aa), and
filters out those overlapping with known CDS annotations from a GFF3 file.

Main class
----------
SixFrameTranslator
    Full pipeline: load genome → six-frame translate → extract ORFs →
    length-filter → overlap-filter → export BED + peptide FASTA.

Typical usage::

    from layer1_sixframe.sixframe import SixFrameTranslator

    st = SixFrameTranslator(
        genome_fasta="genome.fa",
        cds_gff="annotations.gff3",
        out_bed="sorfs.bed",
        out_pep_fasta="sorfs_pep.fasta",
        tmp_dir=None,
        n_jobs=8,
    )
    st.run()
"""

from .sixframe import SixFrameTranslator

__all__ = ["SixFrameTranslator"]
