"""
PeptSesame — Layer4: Functional Prediction & Visualization
==========================================================

Core module providing the :class:`FunctionalAnnotator` which adds
biological context to sORF-encoded peptides after Layer 1 (six-frame),
Layer 2 (scoring), and Layer 3 (classification).

Functional channels
-------------------
1. **GO annotation**       — Homology-based Gene Ontology term transfer.
2. **KEGG mapping**        — BLASTP-based KEGG pathway assignment.
3. **Stress association**  — RNA-seq differential-expression integration.
4. **Receptor pairing**    — LRR-RLK partner prediction for SSP families.
5. **Tissue specificity**  — Multi-tissue expression pattern analysis.
6. **Conservation profile**— Cross-species conservation heatmap.
7. **Visual summary**      — Matplotlib summary plots (histograms, pies, heatmaps).

References
----------
- Tavormina et al. (2015) The Plant Journal
- Olsson et al. (2019) BMC Genomics
- Gholami et al. (2023) Frontiers in Plant Science
"""

from __future__ import annotations

import csv
import itertools
import json
import logging
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ============================================================================
# GO term constants — plant peptide-relevant curated sets
# ============================================================================

#: Broad hierarchy GO terms relevant to small signalling peptides.
BASE_GO_TERMS: Dict[str, Dict[str, Any]] = {
    "GO:0005576": {
        "term": "extracellular region",
        "aspect": "cellular_component",
        "description": "The space external to the plasma membrane — "
                        "canonical location for secreted SSPs.",
    },
    "GO:0005615": {
        "term": "extracellular space",
        "aspect": "cellular_component",
        "description": "Space external to the outermost structure of the cell.",
    },
    "GO:0038023": {
        "term": "signaling receptor binding",
        "aspect": "molecular_function",
        "description": "Binding to a signalling receptor.",
    },
    "GO:0005515": {
        "term": "protein binding",
        "aspect": "molecular_function",
        "description": "Interacting selectively with any protein or peptide.",
    },
    "GO:0007165": {
        "term": "signal transduction",
        "aspect": "biological_process",
        "description": "Cellular process of relaying a signal.",
    },
    "GO:0010469": {
        "term": "regulation of signaling receptor activity",
        "aspect": "biological_process",
        "description": "Modulating signalling receptor activity.",
    },
    "GO:0006952": {
        "term": "defense response",
        "aspect": "biological_process",
        "description": "Reactions triggered by a threat to the organism.",
    },
    "GO:0009733": {
        "term": "response to auxin",
        "aspect": "biological_process",
        "description": "Response to auxin stimulus.",
    },
    "GO:0009751": {
        "term": "response to salicylic acid",
        "aspect": "biological_process",
        "description": "Response to salicylic acid stimulus.",
    },
    "GO:0009753": {
        "term": "response to jasmonic acid",
        "aspect": "biological_process",
        "description": "Response to jasmonic acid stimulus.",
    },
    "GO:0009409": {
        "term": "response to cold",
        "aspect": "biological_process",
        "description": "Response to cold stress.",
    },
    "GO:0009414": {
        "term": "response to water deprivation",
        "aspect": "biological_process",
        "description": "Response to drought / water deprivation.",
    },
    "GO:0050832": {
        "term": "defense response to fungus",
        "aspect": "biological_process",
        "description": "Defence response to fungal infection.",
    },
    "GO:0009651": {
        "term": "response to salt stress",
        "aspect": "biological_process",
        "description": "Response to high salinity.",
    },
    "GO:0006950": {
        "term": "response to stress",
        "aspect": "biological_process",
        "description": "General response to stress stimulus.",
    },
    "GO:0007267": {
        "term": "cell-cell signaling",
        "aspect": "biological_process",
        "description": "Direct intercellular signalling.",
    },
    "GO:0042742": {
        "term": "defense response to bacterium",
        "aspect": "biological_process",
        "description": "Defence response to bacterial infection.",
    },
    "GO:0008544": {
        "term": "epidermis development",
        "aspect": "biological_process",
        "description": "Process of epidermis formation.",
    },
    "GO:0010073": {
        "term": "meristem maintenance",
        "aspect": "biological_process",
        "description": "Maintenance of meristem identity and activity.",
    },
    "GO:0048364": {
        "term": "root development",
        "aspect": "biological_process",
        "description": "Process of root development.",
    },
    "GO:0010223": {
        "term": "secondary shoot formation",
        "aspect": "biological_process",
        "description": "Process of secondary shoot / lateral branch formation.",
    },
}

#: Plant-specific GO terms enriched in peptide datasets.
PLANT_GO_TERMS: Dict[str, Dict[str, Any]] = {
    # Cell-wall-related
    "GO:0042546": {
        "term": "cell wall biogenesis",
        "aspect": "biological_process",
        "description": "Synthesis and assembly of cell wall components.",
    },
    "GO:0009834": {
        "term": "plant-type secondary cell wall biogenesis",
        "aspect": "biological_process",
        "description": "Secondary cell wall formation in plant cells.",
    },
    "GO:0009826": {
        "term": "unidimensional cell growth",
        "aspect": "biological_process",
        "description": "Cell expansion in one dimension, e.g. pollen tube.",
    },
    "GO:0009860": {
        "term": "pollen tube growth",
        "aspect": "biological_process",
        "description": "Growth of the pollen tube.",
    },
    "GO:0010200": {
        "term": "response to chitin",
        "aspect": "biological_process",
        "description": "Response to chitin stimulus.",
    },
    "GO:0009696": {
        "term": "salicylic acid metabolic process",
        "aspect": "biological_process",
        "description": "SA metabolism.",
    },
    "GO:0009744": {
        "term": "response to sucrose",
        "aspect": "biological_process",
        "description": "Response to sucrose stimulus.",
    },
    "GO:0010025": {
        "term": "wax biosynthetic process",
        "aspect": "biological_process",
        "description": "Cuticular wax biosynthesis.",
    },
    "GO:0010075": {
        "term": "regulation of meristem growth",
        "aspect": "biological_process",
        "description": "M, WUS, CLV3-related meristem regulation.",
    },
    "GO:0010154": {
        "term": "fruit development",
        "aspect": "biological_process",
        "description": "Process of fruit formation.",
    },
    "GO:0048440": {
        "term": "carpel development",
        "aspect": "biological_process",
        "description": "Carpel / gynoecium development.",
    },
    "GO:0010582": {
        "term": "floral organ abscission",
        "aspect": "biological_process",
        "description": "Process of floral organ abscission (IDA-related).",
    },
    "GO:0010053": {
        "term": "lateral root formation",
        "aspect": "biological_process",
        "description": "Lateral root development process.",
    },
    "GO:0048441": {
        "term": "stamen development",
        "aspect": "biological_process",
        "description": "Stamen / anther development.",
    },
    "GO:0030154": {
        "term": "cell differentiation",
        "aspect": "biological_process",
        "description": "Cell differentiation process.",
    },
    "GO:0046777": {
        "term": "protein autophosphorylation",
        "aspect": "molecular_function",
        "description": "Autophosphorylation; common in LRR-RLK signalling.",
    },
}

#: Combined GO database for lookups
_ALL_GO_TERMS: Dict[str, Dict[str, Any]] = {}
_ALL_GO_TERMS.update(BASE_GO_TERMS)
_ALL_GO_TERMS.update(PLANT_GO_TERMS)


# ============================================================================
# KEGG pathway constants — plant pathway mappings
# ============================================================================

#: Curated KEGG pathway IDs relevant to plant peptide signalling.
#: Each entry maps pathway_id → {name, category, description}.
KEGG_PLANT_PATHWAYS: Dict[str, Dict[str, str]] = {
    "ath04016": {
        "name": "MAPK signaling pathway — plant",
        "category": "Signal Transduction",
        "description": "Mitogen-activated protein kinase cascade in plants; "
                        "peptide ligands (e.g. PEP1, systemin) activate "
                        "MAPK cascades via receptor kinases.",
    },
    "ath04075": {
        "name": "Plant hormone signal transduction",
        "category": "Signal Transduction",
        "description": "Hormone signalling including auxin, BR, JA, SA; "
                        "some SSPs crosstalk with hormonal pathways.",
    },
    "ath04626": {
        "name": "Plant-pathogen interaction",
        "category": "Environmental Adaptation",
        "description": "Plant innate immunity; RALF, PEP, IDA peptides "
                        "function in defence signalling.",
    },
    "ath04625": {
        "name": "C-type lectin receptor signaling pathway",
        "category": "Signal Transduction",
        "description": "Lectin-RLK mediated signalling in immunity.",
    },
    "ath04370": {
        "name": "VEGF signaling pathway",
        "category": "Signal Transduction",
        "description": "Vascular endothelial growth factor-like signalling; "
                        "related to PSK family.",
    },
    "ath04350": {
        "name": "TGF-beta signaling pathway",
        "category": "Signal Transduction",
        "description": "TGF-beta receptor signalling; shares features "
                        "with CLE/CLV signalling.",
    },
    "ath04141": {
        "name": "Endocytosis",
        "category": "Transport and Catabolism",
        "description": "Endocytic pathway; receptor-mediated endocytosis "
                        "of peptide-bound RLKs.",
    },
    "ath04144": {
        "name": "Autophagy — other",
        "category": "Cellular Processes",
        "description": "Autophagy pathway; stress-induced autophagy "
                        "interfaces with peptide signalling.",
    },
    "ath04712": {
        "name": "Circadian rhythm — plant",
        "category": "Environmental Adaptation",
        "description": "Plant circadian clock; some peptides show "
                        "diurnal expression patterns.",
    },
    "ath03013": {
        "name": "Nucleocytoplasmic transport",
        "category": "Genetic Information Processing",
        "description": "Nuclear transport; relevant for miPEPs acting "
                        "in the nucleus.",
    },
    "ath04151": {
        "name": "PI3K-Akt signaling pathway",
        "category": "Signal Transduction",
        "description": "PI3K/Akt/mTOR signalling; integrates nutrient "
                        "and growth signals.",
    },
    "ath04622": {
        "name": "RIG-I-like receptor signaling pathway",
        "category": "Signal Transduction",
        "description": "RLR-mediated innate immunity; shares signalling "
                        "components with plant PTI.",
    },
    "ath04810": {
        "name": "Regulation of actin cytoskeleton",
        "category": "Cellular Processes",
        "description": "Actin dynamics; RALF signalling modulates "
                        "cytoskeleton remodelling.",
    },
    "ath00940": {
        "name": "Phenylpropanoid biosynthesis",
        "category": "Metabolism",
        "description": "Secondary metabolite biosynthesis; linked to "
                        "cell-wall peptide responses.",
    },
    "ath00941": {
        "name": "Flavonoid biosynthesis",
        "category": "Metabolism",
        "description": "Flavonoid pathway; stress-induced flavonoids "
                        "co-regulate with peptide signals.",
    },
}


# ============================================================================
# SSP–Receptor pairings (curated from literature)
# ============================================================================

#: Known SSP–LRR-RLK receptor pairings.
#: SSP families (keys) map to lists of receptor candidates.
SSP_RECEPTOR_PAIRS: Dict[str, List[Dict[str, Any]]] = {
    "CLE": [
        {
            "receptor": "CLV1",
            "full_name": "CLAVATA 1",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "12628464",
            "confidence": "high",
        },
        {
            "receptor": "CLV2",
            "full_name": "CLAVATA 2",
            "type": "LRR-RLP",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "12628464",
            "confidence": "high",
        },
        {
            "receptor": "BAM1",
            "full_name": "BARELY ANY MERISTEM 1",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "16543447",
            "confidence": "high",
        },
        {
            "receptor": "BAM2",
            "full_name": "BARELY ANY MERISTEM 2",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "16543447",
            "confidence": "medium",
        },
        {
            "receptor": "RPK2",
            "full_name": "RECEPTOR-LIKE PROTEIN KINASE 2",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "18662614",
            "confidence": "medium",
        },
        {
            "receptor": "SOL2",
            "full_name": "SUPPRESSOR OF LLP 2",
            "type": "LRR-RLP",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "18662614",
            "confidence": "medium",
        },
    ],
    "RALF": [
        {
            "receptor": "FERONIA",
            "full_name": "FERONIA",
            "type": "Catharanthus roseus RLK1-like (CrRLK1L)",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "25757470",
            "confidence": "high",
        },
        {
            "receptor": "ANJEA",
            "full_name": "ANJEA",
            "type": "CrRLK1L",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "30726217",
            "confidence": "medium",
        },
        {
            "receptor": "HERCULES1",
            "full_name": "HERCULES 1",
            "type": "CrRLK1L",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "30726217",
            "confidence": "medium",
        },
        {
            "receptor": "BUPS1",
            "full_name": "BUDDHA'S PAPER SEAL 1",
            "type": "CrRLK1L",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "30726217",
            "confidence": "high",
        },
        {
            "receptor": "BUPS2",
            "full_name": "BUDDHA'S PAPER SEAL 2",
            "type": "CrRLK1L",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "30726217",
            "confidence": "high",
        },
    ],
    "CEP": [
        {
            "receptor": "CEPR1",
            "full_name": "CEP RECEPTOR 1",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "25149148",
            "confidence": "high",
        },
        {
            "receptor": "CEPR2",
            "full_name": "CEP RECEPTOR 2",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "25149148",
            "confidence": "high",
        },
    ],
    "PSK": [
        {
            "receptor": "PSKR1",
            "full_name": "PHYTOSULFOKINE RECEPTOR 1",
            "type": "LRR-RLK X",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "8911680",
            "confidence": "high",
        },
        {
            "receptor": "PSKR2",
            "full_name": "PHYTOSULFOKINE RECEPTOR 2",
            "type": "LRR-RLK X",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "8911680",
            "confidence": "high",
        },
    ],
    "PSY1": [
        {
            "receptor": "PSY1R",
            "full_name": "PSY1 RECEPTOR",
            "type": "LRR-RLK X",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "18256640",
            "confidence": "high",
        },
    ],
    "IDA": [
        {
            "receptor": "HAE",
            "full_name": "HAESA",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "18248859",
            "confidence": "high",
        },
        {
            "receptor": "HSL2",
            "full_name": "HAESA-LIKE 2",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "18248859",
            "confidence": "high",
        },
    ],
    "EPFL": [
        {
            "receptor": "ERECTA",
            "full_name": "ERECTA",
            "type": "LRR-RLK XIII",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "26445712",
            "confidence": "high",
        },
        {
            "receptor": "ERL1",
            "full_name": "ERECTA-LIKE 1",
            "type": "LRR-RLK XIII",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "26445712",
            "confidence": "high",
        },
        {
            "receptor": "ERL2",
            "full_name": "ERECTA-LIKE 2",
            "type": "LRR-RLK XIII",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "26445712",
            "confidence": "high",
        },
        {
            "receptor": "TMM",
            "full_name": "TOO MANY MOUTHS",
            "type": "LRR-RLP",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "26445712",
            "confidence": "high",
        },
    ],
    "RGF": [
        {
            "receptor": "RGFR1",
            "full_name": "RGF RECEPTOR 1",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "27040521",
            "confidence": "high",
        },
        {
            "receptor": "RGFR2",
            "full_name": "RGF RECEPTOR 2",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "27040521",
            "confidence": "high",
        },
        {
            "receptor": "RGFR3",
            "full_name": "RGF RECEPTOR 3",
            "type": "LRR-RLK XI",
            "species": "Arabidopsis thaliana",
            "pubmed_id": "27040521",
            "confidence": "high",
        },
    ],
}

#: Default organism code for KEGG queries.
DEFAULT_KEGG_ORG = "ath"  # Arabidopsis thaliana


# ============================================================================
# Data structures
# ============================================================================


@dataclass
class GOAnnotation:
    """A single Gene Ontology term assignment for an sORF-encoded peptide.

    Parameters
    ----------
    go_id : str
        GO identifier (e.g. ``"GO:0005576"``).
    term : str
        Human-readable term name (e.g. ``"extracellular region"``).
    aspect : str
        GO aspect: ``"biological_process"``, ``"molecular_function"``,
        or ``"cellular_component"``.
    evidence : str
        Evidence code (e.g. ``"IEA"``, ``"ISS"``, ``"RCA"``).
    confidence : float
        Confidence score for this assignment in [0, 1].
    source : str
        Source method: ``"homology"``, ``"interproscan"``, ``"blast"``.
    description : str
        Full description of the GO term.
    """

    go_id: str
    term: str
    aspect: str
    evidence: str = "IEA"
    confidence: float = 0.5
    source: str = "homology"
    description: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "GOAnnotation":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class KEGGResult:
    """A single KEGG pathway assignment for an sORF-encoded peptide.

    Parameters
    ----------
    pathway_id : str
        KEGG pathway identifier (e.g. ``"ath04016"``).
    name : str
        Human-readable pathway name.
    category : str
        Pathway category (e.g. ``"Signal Transduction"``).
    confidence : float
        Confidence score in [0, 1].
    description : str
        Pathway description.
    best_hit : str
        Best BLAST hit gene name / locus used for assignment.
    """

    pathway_id: str
    name: str
    category: str
    confidence: float = 0.5
    description: str = ""
    best_hit: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "KEGGResult":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class ReceptorPairingResult:
    """Predicted receptor–ligand pair for an SSP.

    Parameters
    ----------
    ssp_family : str
        SSP family (``"CLE"``, ``"RALF"``, …).
    peptide_id : str
        sORF identifier.
    receptor : str
        Predicted receptor gene name (e.g. ``"CLV1"``).
    full_name : str
        Full receptor name.
    receptor_type : str
        Receptor type (e.g. ``"LRR-RLK XI"``).
    confidence : str
        Confidence string (``"high"``, ``"medium"``, ``"low"``).
    pubmed_id : str
        Supporting literature reference.
    species : str
        Species where the pair was originally described.
    score : float
        Numeric score in [0, 1].
    """

    ssp_family: str
    peptide_id: str
    receptor: str
    full_name: str = ""
    receptor_type: str = ""
    confidence: str = "medium"
    pubmed_id: str = ""
    species: str = "Arabidopsis thaliana"
    score: float = 0.5

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ReceptorPairingResult":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class StressAssociationResult:
    """Result of associating an sORF with a stress condition.

    Parameters
    ----------
    condition : str
        Stress condition name (e.g. ``"drought"``, ``"cold"``).
    log2fc : float
        Log2 fold-change in expression under stress.
    p_value : float
        Adjusted p-value for differential expression.
    status : str
        ``"up"``, ``"down"``, or ``"ns"`` (not significant).
    confidence : float
        Confidence in [0, 1].
    source : str
        Source evidence: ``"rnaseq"``, ``"curated"``, etc.
    """

    condition: str
    log2fc: float = 0.0
    p_value: float = 1.0
    status: str = "ns"
    confidence: float = 0.0
    source: str = "curated"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "StressAssociationResult":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ============================================================================
# Helper utilities
# ============================================================================


def _validate_sequence(seq: str) -> str:
    """Strip whitespace, uppercase, and check for invalid characters.

    Parameters
    ----------
    seq : str
        Amino acid sequence.

    Returns
    -------
    str
        Cleaned, uppercased sequence.

    Raises
    ------
    ValueError
        If the sequence contains invalid amino acid characters.
    """
    VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.strip().upper()
    invalid = set(seq) - VALID_AA
    if invalid:
        raise ValueError(
            f"Invalid amino acid characters in sequence: {''.join(sorted(invalid))}"
        )
    return seq


def _compute_kmer_diversity(seq: str, k: int = 2) -> float:
    """Compute Shannon entropy of k-mer distribution.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    k : int
        k-mer length (default 2).

    Returns
    -------
    float
        Shannon entropy — higher values indicate more diverse composition.
    """
    seq = _validate_sequence(seq)
    if len(seq) < k:
        return 0.0
    total = len(seq) - k + 1
    counter: Counter[str] = Counter()
    for i in range(total):
        counter[seq[i : i + k]] += 1
    entropy = -sum((c / total) * math.log2(c / total) for c in counter.values())
    return entropy


# ============================================================================
# Core annotation class
# ============================================================================


class FunctionalAnnotator:
    """Functional prediction and visualisation for sORF-encoded peptides.

    Adds biological context to candidate small peptides by performing:

    - Homology-based Gene Ontology (GO) term transfer.
    - KEGG pathway assignment via BLASTP.
    - Stress-response association (RNA-seq DE data integration).
    - Receptor kinase pairing for SSP families.
    - Tissue-specificity analysis.
    - Cross-species conservation profiling.
    - Summary visualisation (histograms, pie charts, heatmaps).

    Parameters
    ----------
    blast_db : str or Path, optional
        Path to BLASTP database for homology-based annotation.
    kegg_org : str
        KEGG organism code for pathway mapping (default ``"ath"``).
    expression_data : str or Path, optional
        Path to RNA-seq expression matrix (TSV: gene_id / transcript_id x sample).
    output_dir : str or Path
        Directory for output plots and reports (default: ``"./layer4_output"``).
    n_jobs : int
        Number of parallel worker processes (default: 1).
    """

    def __init__(
        self,
        blast_db: Optional[str | Path] = None,
        kegg_org: str = DEFAULT_KEGG_ORG,
        expression_data: Optional[str | Path] = None,
        output_dir: str | Path = "./layer4_output",
        n_jobs: int = 1,
    ) -> None:
        self.blast_db = Path(blast_db) if blast_db else None
        self.kegg_org = kegg_org
        self.expression_data = Path(expression_data) if expression_data else None
        self.output_dir = Path(output_dir)
        self.n_jobs = n_jobs

        # Internal caches
        self._blast_results: Dict[str, List[Dict[str, Any]]] = {}
        self._expression_matrix: Optional[Dict[str, Dict[str, float]]] = None
        self._kegg_gene_map: Dict[str, List[str]] = {}  # blast_hit → pathway_ids

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            "FunctionalAnnotator initialized (blast_db=%s, kegg_org=%s, "
            "output=%s, n_jobs=%d)",
            self.blast_db,
            self.kegg_org,
            self.output_dir,
            self.n_jobs,
        )

    # ------------------------------------------------------------------
    # 1. GO annotation — homology-based GO term transfer
    # ------------------------------------------------------------------

    def go_annotation(
        self,
        sequence: str,
        seq_id: str = "unknown",
        blast_hits: Optional[List[Dict[str, Any]]] = None,
        similarity_threshold: float = 30.0,
    ) -> List[GOAnnotation]:
        """Annotate a peptide sequence with Gene Ontology terms.

        Uses two strategies:

        1. **Homology-based transfer** — if BLASTP results are provided
           or the BLAST DB is configured, GO terms are transferred from
           the best-matching annotated protein.
        2. **Composition-based inference** — sequence properties
           (cysteine richness, signal peptide, dipeptide entropy) are
           used to infer likely GO aspects when no homology is available.

        Parameters
        ----------
        sequence : str
            Amino acid sequence to annotate.
        seq_id : str
            sORF identifier (for logging).
        blast_hits : list of dict, optional
            Pre-computed BLASTP hits. Each dict must contain 'pident',
            'evalue', 'bitscore', 'qcovs', and optionally 'subject_id'.
        similarity_threshold : float
            Minimum percent identity to transfer GO terms (default 30.0).

        Returns
        -------
        list of GOAnnotation
            Annotated GO terms, sorted by confidence descending.
        """
        seq = _validate_sequence(sequence)
        annotations: List[GOAnnotation] = []
        used_homology = False

        # --- Strategy 1: Homology-based transfer from BLAST hits ---
        if blast_hits is not None and len(blast_hits) > 0:
            transferred = self._transfer_go_from_blast(
                blast_hits, similarity_threshold
            )
            annotations.extend(transferred)
            used_homology = len(transferred) > 0
            logger.debug(
                "%s: transferred %d GO terms from BLAST hits",
                seq_id,
                len(transferred),
            )

        if not used_homology and self.blast_db is not None:
            # Attempt live BLAST (stub — real integration requires subprocess)
            logger.info(
                "%s: BLAST DB available (%s) — not yet invoked inline",
                seq_id,
                self.blast_db,
            )

        # --- Strategy 2: Sequence-property-based inference ---
        inferred = self._infer_go_from_sequence(seq)
        annotations.extend(inferred)

        # --- Deduplicate by go_id, keep highest confidence ---
        annotations = self._deduplicate_go(annotations)

        # Sort by confidence descending
        annotations.sort(key=lambda a: a.confidence, reverse=True)

        return annotations

    def _transfer_go_from_blast(
        self,
        blast_hits: List[Dict[str, Any]],
        similarity_threshold: float = 30.0,
    ) -> List[GOAnnotation]:
        """Transfer GO terms from the best BLAST hit above threshold.

        Parameters
        ----------
        blast_hits : list of dict
            BLASTP hit records.
        similarity_threshold : float
            Minimum percent identity for transfer.

        Returns
        -------
        list of GOAnnotation
            Transferred GO terms.
        """
        annotations: List[GOAnnotation] = []

        # Filter hits above threshold
        good_hits = [
            h
            for h in blast_hits
            if h.get("pident", 0) >= similarity_threshold
            and h.get("evalue", 1) <= 1e-3
        ]
        if not good_hits:
            return annotations

        # Use the best hit
        best = max(good_hits, key=lambda h: h.get("bitscore", 0))
        pident = best.get("pident", 0)
        confidence = min(pident / 100.0, 1.0)

        # Transfer a generic set of GO terms based on identity level
        # In a real implementation, this would query the subject's GO annotation DB.
        # Here we use confidence-weighted defaults.
        if pident >= 70.0:
            # High identity — assign signalling and extracellular terms
            for go_id in ["GO:0005576", "GO:0007165", "GO:0038023"]:
                entry = _ALL_GO_TERMS.get(go_id)
                if entry:
                    annotations.append(
                        GOAnnotation(
                            go_id=go_id,
                            term=entry["term"],
                            aspect=entry["aspect"],
                            evidence="ISS",
                            confidence=round(confidence, 4),
                            source="homology",
                            description=entry.get("description", ""),
                        )
                    )
        elif pident >= 50.0:
            # Moderate identity — assign broader terms
            for go_id in ["GO:0005576", "GO:0007165"]:
                entry = _ALL_GO_TERMS.get(go_id)
                if entry:
                    annotations.append(
                        GOAnnotation(
                            go_id=go_id,
                            term=entry["term"],
                            aspect=entry["aspect"],
                            evidence="ISS",
                            confidence=round(confidence * 0.8, 4),
                            source="homology",
                            description=entry.get("description", ""),
                        )
                    )
        else:
            # Low identity — minimal transfer
            entry = _ALL_GO_TERMS.get("GO:0005576")
            if entry:
                annotations.append(
                    GOAnnotation(
                        go_id="GO:0005576",
                        term=entry["term"],
                        aspect=entry["aspect"],
                        evidence="IEA",
                        confidence=round(confidence * 0.5, 4),
                        source="homology",
                        description=entry.get("description", ""),
                    )
                )

        return annotations

    @staticmethod
    def _infer_go_from_sequence(seq: str) -> List[GOAnnotation]:
        """Infer GO terms from sequence properties when no homology exists.

        Uses simple sequence heuristics:

        - **Cysteine-rich** (>= 4 Cys) → defence / stress response
        - **Hydrophobic N-terminus** → extracellular (signal peptide)
        - **High dipeptide diversity** → signalling / binding
        - **Basic, short** → antimicrobial / defence
        - **Low complexity** → structural / cell wall

        Parameters
        ----------
        seq : str
            Validated amino acid sequence (uppercase).

        Returns
        -------
        list of GOAnnotation
            Inferred GO terms.
        """
        annotations: List[GOAnnotation] = []
        n = len(seq)
        if n == 0:
            return annotations

        # Cysteine count
        cys_count = seq.count("C")
        cys_ratio = cys_count / n if n > 0 else 0.0

        # Hydrophobic N-terminus (signal-peptide-like)
        n_term_hydrophobic = 0
        if n >= 15:
            first15 = seq[:15]
            hydrophobic_aa = set("AILMFWVG")
            n_term_hydrophobic = sum(1 for aa in first15 if aa in hydrophobic_aa)

        # Dipeptide diversity
        diversity = _compute_kmer_diversity(seq, k=2)
        # Max possible entropy for 20 aa dipeptides ≈ log2(400) ≈ 8.64
        norm_diversity = min(diversity / 8.64, 1.0)

        # Net charge (basic residues)
        basic = seq.count("K") + seq.count("R")
        acidic = seq.count("D") + seq.count("E")
        net_charge = basic - acidic

        # --- Heuristic rules ---

        # Rule 1: Cys-rich → defence response (typically CRPs)
        if cys_count >= 4 or cys_ratio >= 0.04:
            # High cysteine → small cysteine-rich protein / defence
            entry = _ALL_GO_TERMS.get("GO:0006952")
            if entry:
                conf = min(0.3 + cys_count * 0.05, 0.85)
                annotations.append(
                    GOAnnotation(
                        go_id="GO:0006952",
                        term=entry["term"],
                        aspect=entry["aspect"],
                        evidence="RCA",
                        confidence=round(conf, 4),
                        source="composition",
                        description=entry.get("description", ""),
                    )
                )

        # Rule 2: Signal-peptide-like N-terminus → extracellular
        if n_term_hydrophobic >= 6 and n >= 20:
            entry = _ALL_GO_TERMS.get("GO:0005576")
            if entry:
                conf = min(0.3 + n_term_hydrophobic * 0.04, 0.80)
                annotations.append(
                    GOAnnotation(
                        go_id="GO:0005576",
                        term=entry["term"],
                        aspect=entry["aspect"],
                        evidence="RCA",
                        confidence=round(conf, 4),
                        source="composition",
                        description=entry.get("description", ""),
                    )
                )

        # Rule 3: High diversity → signalling / binding
        if norm_diversity > 0.5:
            entry = _ALL_GO_TERMS.get("GO:0005515")
            if entry:
                conf = min(0.3 + norm_diversity * 0.3, 0.70)
                annotations.append(
                    GOAnnotation(
                        go_id="GO:0005515",
                        term=entry["term"],
                        aspect=entry["aspect"],
                        evidence="RCA",
                        confidence=round(conf, 4),
                        source="composition",
                        description=entry.get("description", ""),
                    )
                )

        # Rule 4: Basic, short → antimicrobial / defence
        if net_charge >= 2 and n <= 50:
            for go_id in ["GO:0050832", "GO:0042742"]:
                entry = _ALL_GO_TERMS.get(go_id)
                if entry:
                    conf = min(0.3 + net_charge * 0.05, 0.70)
                    annotations.append(
                        GOAnnotation(
                            go_id=go_id,
                            term=entry["term"],
                            aspect=entry["aspect"],
                            evidence="RCA",
                            confidence=round(conf, 4),
                            source="composition",
                            description=entry.get("description", ""),
                        )
                    )

        # Rule 5: Very small (< 20 aa) → possible signalling role
        if n < 20:
            entry = _ALL_GO_TERMS.get("GO:0007267")
            if entry:
                annotations.append(
                    GOAnnotation(
                        go_id="GO:0007267",
                        term=entry["term"],
                        aspect=entry["aspect"],
                        evidence="RCA",
                        confidence=0.3,
                        source="composition",
                        description=entry.get("description", ""),
                    )
                )

        return annotations

    @staticmethod
    def _deduplicate_go(
        annotations: List[GOAnnotation],
    ) -> List[GOAnnotation]:
        """Deduplicate GO annotations by go_id, keeping highest confidence.

        Parameters
        ----------
        annotations : list of GOAnnotation
            Input annotations (possibly with duplicates).

        Returns
        -------
        list of GOAnnotation
            Deduplicated annotations.
        """
        seen: Dict[str, GOAnnotation] = {}
        for ann in annotations:
            if ann.go_id not in seen or ann.confidence > seen[ann.go_id].confidence:
                seen[ann.go_id] = ann
        return list(seen.values())

    # ------------------------------------------------------------------
    # 2. KEGG pathway mapping
    # ------------------------------------------------------------------

    def kegg_mapping(
        self,
        sequence: str,
        seq_id: str = "unknown",
        blast_hits: Optional[List[Dict[str, Any]]] = None,
    ) -> List[KEGGResult]:
        """Assign KEGG pathways to a peptide sequence.

        Uses BLASTP homology to known KEGG-annotated proteins.
        When no BLAST hit is available, uses sequence-property heuristics
        to suggest likely pathway categories.

        Parameters
        ----------
        sequence : str
            Amino acid sequence.
        seq_id : str
            sORF identifier.
        blast_hits : list of dict, optional
            Pre-computed BLASTP hits (each dict should contain 'pident',
            'evalue', 'bitscore', and optionally 'subject_id' and 'kegg_gene').

        Returns
        -------
        list of KEGGResult
            Assigned KEGG pathways, sorted by confidence descending.
        """
        seq = _validate_sequence(sequence)
        results: List[KEGGResult] = []
        used_homology = False

        # --- Strategy 1: Homology-based via BLAST hits ---
        if blast_hits is not None and len(blast_hits) > 0:
            mapped = self._map_kegg_from_blast(blast_hits)
            results.extend(mapped)
            used_homology = len(mapped) > 0

        # --- Strategy 2: Composition-based heuristic inference ---
        if not used_homology:
            inferred = self._infer_kegg_from_sequence(seq)
            results.extend(inferred)

        # Deduplicate by pathway_id, keep highest confidence
        results = self._deduplicate_kegg(results)
        results.sort(key=lambda r: r.confidence, reverse=True)

        return results

    @staticmethod
    def _map_kegg_from_blast(
        blast_hits: List[Dict[str, Any]],
    ) -> List[KEGGResult]:
        """Map KEGG pathways from BLAST hits.

        Parameters
        ----------
        blast_hits : list of dict
            BLASTP hit records.

        Returns
        -------
        list of KEGGResult
            Mapped KEGG pathways.
        """
        results: List[KEGGResult] = []

        # Filter hits
        good = [
            h
            for h in blast_hits
            if h.get("evalue", 1) <= 1e-3 and h.get("pident", 0) >= 30.0
        ]
        if not good:
            return results

        best = max(good, key=lambda h: h.get("bitscore", 0))
        pident = best.get("pident", 0)
        subject_id = best.get("subject_id", best.get("kegg_gene", "unknown"))
        confidence = min(pident / 100.0, 1.0)

        # Map known KEGG-annotated genes to pathways
        # In a real pipeline, this would query the KEGG API or a local DB.
        # Here we use the curated pathway set with the hit gene as context.
        default_pathways = [
            "ath04016",  # MAPK signalling
            "ath04626",  # Plant-pathogen interaction
        ]
        for pid in default_pathways:
            entry = KEGG_PLANT_PATHWAYS.get(pid)
            if entry:
                results.append(
                    KEGGResult(
                        pathway_id=pid,
                        name=entry["name"],
                        category=entry["category"],
                        confidence=round(confidence * 0.8, 4),
                        description=entry["description"],
                        best_hit=subject_id,
                    )
                )

        # If identity is very high, also add hormone signalling
        if pident >= 60.0:
            entry = KEGG_PLANT_PATHWAYS.get("ath04075")
            if entry:
                results.append(
                    KEGGResult(
                        pathway_id="ath04075",
                        name=entry["name"],
                        category=entry["category"],
                        confidence=round(confidence * 0.9, 4),
                        description=entry["description"],
                        best_hit=subject_id,
                    )
                )

        return results

    @staticmethod
    def _infer_kegg_from_sequence(seq: str) -> List[KEGGResult]:
        """Infer likely KEGG pathways from sequence properties.

        Parameters
        ----------
        seq : str
            Validated amino acid sequence.

        Returns
        -------
        list of KEGGResult
            Inferred pathway assignments.
        """
        results: List[KEGGResult] = []
        n = len(seq)
        if n == 0:
            return results

        # Cys content
        cys_ratio = seq.count("C") / n
        # Signal peptide heuristic
        n_term_hydro = sum(1 for aa in seq[:15] if aa in "AILMFWVG") if n >= 15 else 0

        # Rule 1: Cys-rich → plant-pathogen interaction / defence
        if cys_ratio >= 0.04 or seq.count("C") >= 4:
            for pid in ["ath04626", "ath04016"]:
                entry = KEGG_PLANT_PATHWAYS.get(pid)
                if entry:
                    conf = min(0.3 + seq.count("C") * 0.05, 0.7)
                    results.append(
                        KEGGResult(
                            pathway_id=pid,
                            name=entry["name"],
                            category=entry["category"],
                            confidence=round(conf, 4),
                            description=entry["description"],
                            best_hit="inferred",
                        )
                    )

        # Rule 2: Putative secreted → endocytosis / signalling
        if n_term_hydro >= 6 and n >= 20:
            for pid in ["ath04141", "ath04075"]:
                entry = KEGG_PLANT_PATHWAYS.get(pid)
                if entry:
                    results.append(
                        KEGGResult(
                            pathway_id=pid,
                            name=entry["name"],
                            category=entry["category"],
                            confidence=0.3,
                            description=entry["description"],
                            best_hit="inferred",
                        )
                    )

        # Rule 3: Very small (< 15 aa) → possible hormone-like signalling
        if n < 15:
            entry = KEGG_PLANT_PATHWAYS.get("ath04075")
            if entry:
                results.append(
                    KEGGResult(
                        pathway_id="ath04075",
                        name=entry["name"],
                        category=entry["category"],
                        confidence=0.25,
                        description=entry["description"],
                        best_hit="inferred",
                    )
                )

        return results

    @staticmethod
    def _deduplicate_kegg(
        results: List[KEGGResult],
    ) -> List[KEGGResult]:
        """Deduplicate KEGG results by pathway_id, keeping highest confidence.

        Parameters
        ----------
        results : list of KEGGResult
            Input results.

        Returns
        -------
        list of KEGGResult
            Deduplicated results.
        """
        seen: Dict[str, KEGGResult] = {}
        for r in results:
            if r.pathway_id not in seen or r.confidence > seen[r.pathway_id].confidence:
                seen[r.pathway_id] = r
        return list(seen.values())

    # ------------------------------------------------------------------
    # 3. Stress association — RNA-seq DE data integration
    # ------------------------------------------------------------------

    def stress_association(
        self,
        seq_id: str,
        expression_profile: Optional[Dict[str, float]] = None,
        de_results: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> List[StressAssociationResult]:
        """Associate an sORF with stress conditions.

        Operates in three modes:

        1. **DE results provided** — directly use pre-computed differential
           expression statistics.
        2. **Expression profile provided** — compare expression profile
           against known stress marker signatures.
        3. **No data** — return empty list with a placeholder note.

        Parameters
        ----------
        seq_id : str
            sORF identifier.
        expression_profile : dict, optional
            Expression values (TPM / RPKM) keyed by sample/condition.
        de_results : dict, optional
            Differential expression results keyed by condition. Each value
            is a dict with ``"log2FC"`` and ``"pvalue"`` keys.

        Returns
        -------
        list of StressAssociationResult
            Per-condition associations.
        """
        results: List[StressAssociationResult] = []

        if de_results is not None:
            # Mode 1: Direct DE results
            for condition, stats in de_results.items():
                log2fc = stats.get("log2FC", 0.0)
                pval = stats.get("pvalue", 1.0)
                abs_log2fc = abs(log2fc)

                # Determine status
                if pval <= 0.05 and abs_log2fc >= 1.0:
                    status = "up" if log2fc > 0 else "down"
                elif pval <= 0.05:
                    status = "up" if log2fc > 0 else "down"
                else:
                    status = "ns"

                confidence = 0.0
                if status != "ns":
                    # Confidence from effect size and significance
                    conf_log2fc = min(abs_log2fc / 5.0, 1.0)
                    conf_pval = 1.0 - min(pval, 1.0)
                    confidence = round(0.6 * conf_log2fc + 0.4 * conf_pval, 4)
                    confidence = min(max(confidence, 0.1), 0.95)

                results.append(
                    StressAssociationResult(
                        condition=condition,
                        log2fc=round(log2fc, 4),
                        p_value=round(pval, 6),
                        status=status,
                        confidence=confidence,
                        source="rnaseq",
                    )
                )

        elif expression_profile is not None:
            # Mode 2: Expression pattern-based inference
            # Compare against known stress-marker expression signatures
            stress_signatures = self._get_stress_signatures()
            for condition, sig_genes in stress_signatures.items():
                # Compute correlation between profile and signature
                # (simplified: check if any signature gene is present)
                score = self._profile_stress_score(
                    expression_profile, sig_genes
                )
                if score > 0:
                    results.append(
                        StressAssociationResult(
                            condition=condition,
                            log2fc=0.0,
                            p_value=1.0,
                            status="up" if score > 0.3 else "ns",
                            confidence=round(min(score, 0.8), 4),
                            source="expression_pattern",
                        )
                    )

        else:
            # Mode 3: No data — note in logger
            logger.debug("%s: no expression / DE data — skipping stress association", seq_id)

        return results

    @staticmethod
    def _get_stress_signatures() -> Dict[str, List[str]]:
        """Get known stress-response marker gene signatures.

        Returns
        -------
        dict
            Mapping of stress condition → list of marker gene names.
        """
        return {
            "drought": ["RD29A", "RD29B", "LEA", "DREB2A", "DREB2B"],
            "cold": ["CBF1", "CBF2", "CBF3", "COR15A", "COR47"],
            "salt": ["SOS1", "SOS2", "SOS3", "NHX1", "HKT1"],
            "heat": ["HSP70", "HSP90", "HSP101", "HSFA1", "HSFA2"],
            "pathogen": ["PR1", "PR2", "PR5", "PDF1.2", "WRKY33"],
            "wounding": ["JAZ1", "JAZ2", "JAZ3", "OPR3", "AOS"],
            "oxidative": ["APX1", "APX2", "SOD", "CAT1", "GPX"],
        }

    @staticmethod
    def _profile_stress_score(
        profile: Dict[str, float],
        signature_genes: List[str],
    ) -> float:
        """Compute a stress-association score from expression profile.

        Simple heuristic: checks if the expression profile keys overlap
        with known stress condition names or if the profile values
        resemble a stress response pattern.

        Parameters
        ----------
        profile : dict
            Expression values per sample/condition.
        signature_genes : list of str
            Known stress marker genes for a condition.

        Returns
        -------
        float
            Score in [0, 1] indicating strength of association.
        """
        # Check if any profile sample name matches the condition
        # or contains condition-related keywords
        score = 0.0
        profile_keys = " ".join(profile.keys()).lower()

        # Check signature gene names against profile keys
        for marker in signature_genes:
            if marker.lower() in profile_keys:
                score += 0.15

        # Check stress-related keywords
        stress_keywords = [
            "stress", "drought", "cold", "salt", "heat", "pathogen",
            "wound", "oxidative", "chitin", "flagellin", "elf18",
        ]
        for kw in stress_keywords:
            if kw in profile_keys:
                score += 0.1

        return min(score, 1.0)

    # ------------------------------------------------------------------
    # 4. Receptor pairing — LRR-RLK prediction for SSPs
    # ------------------------------------------------------------------

    def receptor_pairing(
        self,
        ssp_family: str,
        peptide_id: str,
        sequence: Optional[str] = None,
        confidence_boost: float = 0.0,
    ) -> List[ReceptorPairingResult]:
        """Predict receptor kinase partners for an SSP.

        Uses the curated ``SSP_RECEPTOR_PAIRS`` dictionary to look up
        known LRR-RLK / CrRLK1L receptors for the given SSP family.
        Optionally refines the score based on sequence features.

        Parameters
        ----------
        ssp_family : str
            SSP family name (uppercase, e.g. ``"CLE"``, ``"RALF"``).
        peptide_id : str
            sORF identifier.
        sequence : str, optional
            Full peptide sequence (used for score refinement).
        confidence_boost : float
            Additional confidence boost (e.g. from Layer 2 score).

        Returns
        -------
        list of ReceptorPairingResult
            Predicted receptor pairings, sorted by score descending.
        """
        family_upper = ssp_family.upper().strip()

        pairs = SSP_RECEPTOR_PAIRS.get(family_upper, [])
        if not pairs:
            logger.debug(
                "No known receptor pairs for SSP family '%s'", ssp_family
            )
            return []

        results: List[ReceptorPairingResult] = []
        seq_length = len(sequence) if sequence else 0

        for pair in pairs:
            # Base score from literature confidence
            conf_str = pair.get("confidence", "medium")
            if conf_str == "high":
                base_score = 0.80
            elif conf_str == "medium":
                base_score = 0.55
            else:
                base_score = 0.30

            # Refine score with sequence features
            if sequence:
                seq = _validate_sequence(sequence)
                n = len(seq)

                # Cys count bonus: CLE peptides typically have 6 Cys
                cys_count = seq.count("C")
                if family_upper == "CLE" and cys_count >= 4:
                    base_score += 0.10
                elif family_upper == "RALF" and cys_count >= 4:
                    base_score += 0.10

                # Length bonus: match typical active peptide length
                if 40 <= n <= 140:
                    base_score += 0.05
                elif 10 <= n <= 39:
                    base_score += 0.02

            # Apply confidence boost
            final_score = min(base_score + confidence_boost, 0.99)

            results.append(
                ReceptorPairingResult(
                    ssp_family=family_upper,
                    peptide_id=peptide_id,
                    receptor=pair["receptor"],
                    full_name=pair["full_name"],
                    receptor_type=pair["type"],
                    confidence=conf_str,
                    pubmed_id=pair.get("pubmed_id", ""),
                    species=pair.get("species", "Arabidopsis thaliana"),
                    score=round(final_score, 4),
                )
            )

        results.sort(key=lambda r: r.score, reverse=True)
        return results

    # ------------------------------------------------------------------
    # 5. Tissue specificity
    # ------------------------------------------------------------------

    def tissue_specificity(
        self,
        expression_profile: Dict[str, float],
        tissue_groups: Optional[Dict[str, List[str]]] = None,
    ) -> Dict[str, Any]:
        """Compute tissue-specificity metrics from expression data.

        Uses the Tau index (Yanai et al., 2005) to quantify tissue
        specificity, plus per-group enrichments.

        Parameters
        ----------
        expression_profile : dict
            Expression values (TPM) keyed by sample/tissue name.
        tissue_groups : dict, optional
            Mapping of tissue group → list of sample names in that group.
            If not provided, each sample is treated as its own group.

        Returns
        -------
        dict with keys:
            - ``"tau"`` : float — Tau specificity index [0=ubiquitous, 1=specific].
            - ``"max_tissue"`` : str — Tissue with highest expression.
            - ``"max_value"`` : float — Maximum expression value.
            - ``"enrichment"`` : dict — Tissue → enrichment score.
            - ``"profile"`` : dict — Tissue → normalised expression.
        """
        if not expression_profile:
            return {
                "tau": 0.0,
                "max_tissue": "",
                "max_value": 0.0,
                "enrichment": {},
                "profile": {},
            }

        # --- Normalise: max expression = 1.0 ---
        values = np.array(list(expression_profile.values()), dtype=float)
        tissues = list(expression_profile.keys())
        max_val = float(np.max(values)) if values.size > 0 else 0.0
        normalised = values / max_val if max_val > 0 else np.zeros_like(values)

        profile_dict = dict(zip(tissues, normalised.tolist()))

        # --- Tau index ---
        n = len(normalised)
        if n <= 1:
            tau = 0.0
        else:
            tau = float(np.sum(1.0 - normalised) / (n - 1))
        tau = round(max(0.0, min(1.0, tau)), 4)

        # --- Per-group enrichment ---
        if tissue_groups is None:
            tissue_groups = {t: [t] for t in tissues}

        enrichment: Dict[str, float] = {}
        for group, members in tissue_groups.items():
            group_vals = [
                expression_profile.get(m, 0.0) for m in members
            ]
            group_mean = np.mean(group_vals) if group_vals else 0.0
            all_mean = np.mean(values) if values.size > 0 else 1.0
            enrichment[group] = round(
                float(group_mean / all_mean) if all_mean > 0 else 0.0, 4
            )

        # Max tissue
        max_idx = int(np.argmax(values)) if values.size > 0 else 0
        max_tissue = tissues[max_idx] if tissues else ""

        return {
            "tau": tau,
            "max_tissue": max_tissue,
            "max_value": round(float(max_val), 4),
            "enrichment": enrichment,
            "profile": profile_dict,
        }

    # ------------------------------------------------------------------
    # 6. Conservation profiling
    # ------------------------------------------------------------------

    def conservation_profile(
        self,
        seq_id: str,
        query_sequence: str,
        blast_results: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    ) -> Dict[str, Any]:
        """Build a cross-species conservation profile for an sORF.

        When BLAST results against multiple species are provided,
        constructs a conservation matrix (species → identity),
        and computes summary statistics.

        Parameters
        ----------
        seq_id : str
            sORF identifier.
        query_sequence : str
            Query amino acid sequence.
        blast_results : dict, optional
            Mapping of species name → list of BLAST hit dicts.
            Each hit dict should contain ``"pident"``, ``"evalue"``,
            ``"bitscore"``, and ``"qcovs"``.

        Returns
        -------
        dict with keys:
            - ``"matrix"`` : dict — Species → percent identity.
            - ``"best_hit"`` : str — Species with highest identity.
            - ``"max_identity"`` : float — Maximum percent identity.
            - ``"n_species_with_hits"`` : int — Species count with hits.
            - ``"conservation_score"`` : float — Aggregate score [0, 1].
            - ``"is_conserved"`` : bool — Whether the sORF is conserved.
        """
        seq = _validate_sequence(query_sequence)

        if blast_results is None:
            return {
                "matrix": {},
                "best_hit": "",
                "max_identity": 0.0,
                "n_species_with_hits": 0,
                "conservation_score": 0.0,
                "is_conserved": False,
            }

        matrix: Dict[str, float] = {}
        best_species = ""
        max_identity = 0.0
        n_with_hits = 0

        for species, hits in blast_results.items():
            if not hits:
                continue

            # Filter to significant hits
            good = [
                h
                for h in hits
                if h.get("evalue", 1) <= 1e-3 and h.get("pident", 0) >= 20.0
            ]
            if not good:
                continue

            # Best hit per species
            best = max(good, key=lambda h: h.get("bitscore", 0))
            pident = best.get("pident", 0.0)
            matrix[species] = round(float(pident), 2)
            n_with_hits += 1

            if pident > max_identity:
                max_identity = pident
                best_species = species

        # Aggregate conservation score
        if n_with_hits == 0:
            conservation_score = 0.0
        else:
            avg_identity = np.mean(list(matrix.values()))
            # Score combines average identity and species breadth
            breadth = n_with_hits / max(len(matrix) if matrix else 1, 1)
            conservation_score = round(
                float(0.6 * (avg_identity / 100.0) + 0.4 * breadth), 4
            )
            conservation_score = min(conservation_score, 1.0)

        return {
            "matrix": matrix,
            "best_hit": best_species,
            "max_identity": round(float(max_identity), 2),
            "n_species_with_hits": n_with_hits,
            "conservation_score": conservation_score,
            "is_conserved": conservation_score >= 0.4,
        }

    # ------------------------------------------------------------------
    # 7. Summary plots
    # ------------------------------------------------------------------

    def generate_plots(
        self,
        annotations: List[Dict[str, Any]],
        output_prefix: str = "layer4_summary",
        dpi: int = 150,
    ) -> Dict[str, str]:
        """Generate summary visualisation plots from functional annotations.

        Creates up to four plot types:

        1. **GO aspect pie chart** — distribution of GO aspects
           (biological_process, molecular_function, cellular_component).
        2. **KEGG pathway bar chart** — top KEGG pathway categories.
        3. **Conservation heatmap** — cross-species identity matrix.
        4. **SSP family distribution** — bar chart of SSP families.

        Requires ``matplotlib``.

        Parameters
        ----------
        annotations : list of dict
            List of sORF annotation dicts. Each dict should contain at
            minimum ``"seq_id"``, and optionally ``"go_annotations"``,
            ``"kegg_results"``, ``"ssp_family"``, and ``"conservation"``.
        output_prefix : str
            Prefix for output image files (default ``"layer4_summary"``).
        dpi : int
            Figure resolution in DPI (default 150).

        Returns
        -------
        dict
            Mapping of plot name → file path (relative to ``output_dir``).
            Keys: ``"go_pie"``, ``"kegg_bar"``, ``"conservation_heatmap"``,
            ``"ssp_bar"``. Unavailable plots have value ``""``.
        """
        try:
            import matplotlib
            matplotlib.use("Agg")  # non-interactive backend
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            logger.warning(
                "matplotlib not available — skipping plot generation"
            )
            return {
                "go_pie": "",
                "kegg_bar": "",
                "conservation_heatmap": "",
                "ssp_bar": "",
            }

        generated: Dict[str, str] = {
            "go_pie": "",
            "kegg_bar": "",
            "conservation_heatmap": "",
            "ssp_bar": "",
        }

        n = len(annotations)
        if n == 0:
            logger.warning("No annotations provided for plotting")
            return generated

        # ----------------------------------------------------------
        # 7a. GO aspect pie chart
        # ----------------------------------------------------------
        go_counts: Dict[str, int] = {
            "biological_process": 0,
            "molecular_function": 0,
            "cellular_component": 0,
        }
        total_go = 0
        for ann in annotations:
            go_list = ann.get("go_annotations", [])
            if isinstance(go_list, list):
                for go_ann in go_list:
                    aspect = go_ann.get("aspect", "")
                    if aspect in go_counts:
                        go_counts[aspect] += 1
                        total_go += 1

        if total_go > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            labels = [k.replace("_", " ").title() for k, v in go_counts.items() if v > 0]
            sizes = [v for v in go_counts.values() if v > 0]
            colors_go = ["#66c2a5", "#fc8d62", "#8da0cb"]
            wedges, texts, autotexts = ax.pie(
                sizes,
                labels=labels,
                autopct="%1.1f%%",
                startangle=90,
                colors=colors_go[: len(sizes)],
                wedgeprops={"edgecolor": "white", "linewidth": 1.5},
            )
            ax.set_title("GO Aspect Distribution", fontsize=14, fontweight="bold")
            plt.tight_layout()
            out_path = self.output_dir / f"{output_prefix}_go_pie.png"
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            generated["go_pie"] = str(out_path)
            logger.info("Saved GO pie chart to %s", out_path)

        # ----------------------------------------------------------
        # 7b. KEGG pathway bar chart (top categories)
        # ----------------------------------------------------------
        kegg_categories: Counter[str] = Counter()
        for ann in annotations:
            kegg_list = ann.get("kegg_results", [])
            if isinstance(kegg_list, list):
                for kr in kegg_list:
                    cat = kr.get("category", "Unknown")
                    kegg_categories[cat] += 1

        if kegg_categories:
            fig, ax = plt.subplots(figsize=(8, 5))
            cats_sorted = kegg_categories.most_common(8)
            cat_names = [c[0] for c in cats_sorted]
            cat_counts = [c[1] for c in cats_sorted]
            colors_kegg = plt.cm.Set2(
                np.linspace(0, 1, len(cat_names))
            )
            bars = ax.barh(
                range(len(cat_names)),
                cat_counts,
                color=colors_kegg,
                edgecolor="grey",
                linewidth=0.5,
            )
            ax.set_yticks(range(len(cat_names)))
            ax.set_yticklabels(cat_names, fontsize=10)
            ax.invert_yaxis()
            ax.set_xlabel("Number of assignments", fontsize=11)
            ax.set_title("KEGG Pathway Categories", fontsize=14, fontweight="bold")
            for bar, count in zip(bars, cat_counts):
                ax.text(
                    bar.get_width() + 0.3,
                    bar.get_y() + bar.get_height() / 2,
                    str(count),
                    va="center",
                    fontsize=9,
                )
            plt.tight_layout()
            out_path = self.output_dir / f"{output_prefix}_kegg_bar.png"
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            generated["kegg_bar"] = str(out_path)
            logger.info("Saved KEGG bar chart to %s", out_path)

        # ----------------------------------------------------------
        # 7c. Conservation heatmap
        # ----------------------------------------------------------
        # Collect conservation matrices from annotations that have them
        species_set: Set[str] = set()
        conserved_ids: List[str] = []
        conservation_rows: List[List[float]] = []

        for ann in annotations:
            cons = ann.get("conservation", {})
            if isinstance(cons, dict):
                matrix = cons.get("matrix", {})
                if matrix:
                    species_set.update(matrix.keys())
                    conserved_ids.append(ann.get("seq_id", "unknown"))

        if species_set and conserved_ids:
            species_list = sorted(species_set)
            # Build matrix
            heatmap_data = np.zeros((len(conserved_ids), len(species_list)))
            for i, ann in enumerate(annotations):
                cons = ann.get("conservation", {})
                matrix = cons.get("matrix", {})
                for j, sp in enumerate(species_list):
                    heatmap_data[i, j] = matrix.get(sp, 0.0)

            if heatmap_data.size > 0 and np.any(heatmap_data > 0):
                fig, ax = plt.subplots(figsize=(10, max(4, len(conserved_ids) * 0.3)))
                norm = Normalize(vmin=0, vmax=100)
                im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd", norm=norm)
                ax.set_xticks(range(len(species_list)))
                ax.set_xticklabels(species_list, rotation=45, ha="right", fontsize=8)
                ax.set_yticks(range(len(conserved_ids)))
                ax.set_yticklabels(conserved_ids, fontsize=7)
                ax.set_xlabel("Species", fontsize=11)
                ax.set_title("Cross-Species Conservation (% identity)", fontsize=12, fontweight="bold")
                cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label("% identity", fontsize=10)
                plt.tight_layout()
                out_path = self.output_dir / f"{output_prefix}_conservation_heatmap.png"
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
                generated["conservation_heatmap"] = str(out_path)
                logger.info("Saved conservation heatmap to %s", out_path)

        # ----------------------------------------------------------
        # 7d. SSP family distribution bar chart
        # ----------------------------------------------------------
        ssp_families: Counter[str] = Counter()
        for ann in annotations:
            ssp = ann.get("ssp_family")
            if ssp and isinstance(ssp, str):
                ssp_families[ssp] += 1
            # Also check categories
            categories = ann.get("categories", {})
            if isinstance(categories, dict):
                for cat, conf in categories.items():
                    if cat == "SSP":
                        # The specific family might be in 'ssp_family'
                        pass

        if ssp_families:
            fig, ax = plt.subplots(figsize=(8, 5))
            families = sorted(ssp_families.keys())
            counts = [ssp_families[f] for f in families]
            colors_ssp = plt.cm.Set3(np.linspace(0, 1, len(families)))
            bars = ax.bar(families, counts, color=colors_ssp, edgecolor="grey", linewidth=0.5)
            ax.set_xlabel("SSP Family", fontsize=11)
            ax.set_ylabel("Count", fontsize=11)
            ax.set_title("SSP Family Distribution", fontsize=14, fontweight="bold")
            for bar, count in zip(bars, counts):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                )
            plt.tight_layout()
            out_path = self.output_dir / f"{output_prefix}_ssp_distribution.png"
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            generated["ssp_bar"] = str(out_path)
            logger.info("Saved SSP distribution bar chart to %s", out_path)

        return generated

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    def annotate_batch(
        self,
        sorf_records: List[Dict[str, Any]],
        run_go: bool = True,
        run_kegg: bool = True,
        run_stress: bool = False,
        run_receptor: bool = True,
        run_conservation: bool = False,
        run_tissue: bool = False,
        run_plots: bool = True,
    ) -> List[Dict[str, Any]]:
        """Run functional annotation on a batch of sORF records.

        Each record in ``sorf_records`` should be a dict with at least
        ``"seq_id"`` and ``"sequence"``. Additional fields (``"blast_hits"``,
        ``"ssp_family"``, ``"expression_profile"``, etc.) enable the
        corresponding analyses.

        Parameters
        ----------
        sorf_records : list of dict
            Batch of sORF records to annotate.
        run_go : bool
            Run GO annotation (default True).
        run_kegg : bool
            Run KEGG mapping (default True).
        run_stress : bool
            Run stress association (default False — requires expression data).
        run_receptor : bool
            Run receptor pairing (default True).
        run_conservation : bool
            Run conservation profiling (default False — requires BLAST results).
        run_tissue : bool
            Run tissue specificity (default False — requires expression profile).
        run_plots : bool
            Generate summary plots (default True).

        Returns
        -------
        list of dict
            Annotated sORF records with functional data added as new keys.
        """
        logger.info(
            "Annotating batch of %d sORFs (go=%s, kegg=%s, stress=%s, "
            "receptor=%s, conservation=%s, tissue=%s, plots=%s)",
            len(sorf_records),
            run_go,
            run_kegg,
            run_stress,
            run_receptor,
            run_conservation,
            run_tissue,
            run_plots,
        )

        annotated: List[Dict[str, Any]] = []

        for record in sorf_records:
            result = dict(record)  # shallow copy
            seq_id = result.get("seq_id", "unknown")
            sequence = result.get("sequence", "")
            ssp_family = result.get("ssp_family")

            if not sequence:
                logger.warning("%s: no sequence — skipping", seq_id)
                annotated.append(result)
                continue

            # 1. GO annotation
            if run_go:
                blast_hits = result.get("blast_hits")
                go_anns = self.go_annotation(
                    sequence=sequence,
                    seq_id=seq_id,
                    blast_hits=blast_hits,
                )
                result["go_annotations"] = [g.to_dict() for g in go_anns]
            else:
                result["go_annotations"] = []

            # 2. KEGG mapping
            if run_kegg:
                blast_hits = result.get("blast_hits")
                kegg_results = self.kegg_mapping(
                    sequence=sequence,
                    seq_id=seq_id,
                    blast_hits=blast_hits,
                )
                result["kegg_results"] = [k.to_dict() for k in kegg_results]
            else:
                result["kegg_results"] = []

            # 3. Stress association
            if run_stress:
                de_results = result.get("de_results")
                expr_profile = result.get("expression_profile")
                stress = self.stress_association(
                    seq_id=seq_id,
                    expression_profile=expr_profile,
                    de_results=de_results,
                )
                result["stress_associations"] = [s.to_dict() for s in stress]
            else:
                result["stress_associations"] = []

            # 4. Receptor pairing
            if run_receptor and ssp_family:
                confidence_boost = result.get("score", 0.0)
                pairs = self.receptor_pairing(
                    ssp_family=ssp_family,
                    peptide_id=seq_id,
                    sequence=sequence,
                    confidence_boost=confidence_boost,
                )
                result["receptor_pairs"] = [p.to_dict() for p in pairs]
            else:
                result["receptor_pairs"] = []

            # 5. Conservation profiling
            if run_conservation:
                blast_by_species = result.get("blast_by_species")
                cons = self.conservation_profile(
                    seq_id=seq_id,
                    query_sequence=sequence,
                    blast_results=blast_by_species,
                )
                result["conservation"] = cons
            else:
                result["conservation"] = {}

            # 6. Tissue specificity
            if run_tissue:
                expr_profile = result.get("expression_profile", {})
                tissue_groups = result.get("tissue_groups")
                ts = self.tissue_specificity(
                    expression_profile=expr_profile,
                    tissue_groups=tissue_groups,
                )
                result["tissue_specificity"] = ts
            else:
                result["tissue_specificity"] = {}

            annotated.append(result)

        # 7. Summary plots (on the full batch)
        if run_plots:
            plots = self.generate_plots(annotated)
            logger.info("Generated plots: %s", plots)
            # Attach plot paths to each record
            for rec in annotated:
                rec["layer4_plots"] = plots

        logger.info("Batch annotation complete — %d records processed", len(annotated))
        return annotated

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    def export_annotations(
        self,
        annotations: List[Dict[str, Any]],
        format: str = "json",
        output_path: Optional[str | Path] = None,
    ) -> str:
        """Export functional annotations to file.

        Parameters
        ----------
        annotations : list of dict
            Annotated sORF records.
        format : str
            Output format: ``"json"`` or ``"tsv"``.
        output_path : str or Path, optional
            Output file path. Defaults to
            ``{output_dir}/layer4_annotations.{format}``.

        Returns
        -------
        str
            Path to the written output file.
        """
        if output_path is None:
            output_path = self.output_dir / f"layer4_annotations.{format}"

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "json":
            with open(output_path, "w") as f:
                json.dump(annotations, f, indent=2, default=str)
            logger.info("Exported JSON annotations to %s", output_path)

        elif format == "tsv":
            # Flatten key annotation fields into TSV
            fieldnames = [
                "seq_id",
                "ssp_family",
                "go_terms",
                "kegg_pathways",
                "receptor_pairs",
                "conservation_score",
                "stress_conditions",
            ]
            with open(output_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for ann in annotations:
                    row = {
                        "seq_id": ann.get("seq_id", ""),
                        "ssp_family": ann.get("ssp_family", ""),
                        "go_terms": ";".join(
                            g.get("go_id", "")
                            for g in ann.get("go_annotations", [])
                        ),
                        "kegg_pathways": ";".join(
                            k.get("pathway_id", "")
                            for k in ann.get("kegg_results", [])
                        ),
                        "receptor_pairs": ";".join(
                            p.get("receptor", "")
                            for p in ann.get("receptor_pairs", [])
                        ),
                        "conservation_score": str(
                            ann.get("conservation", {}).get("conservation_score", "")
                        ),
                        "stress_conditions": ";".join(
                            s.get("condition", "")
                            for s in ann.get("stress_associations", [])
                            if s.get("status") != "ns"
                        ),
                    }
                    writer.writerow(row)
            logger.info("Exported TSV annotations to %s", output_path)

        else:
            raise ValueError(f"Unsupported export format: {format}")

        return str(output_path)
