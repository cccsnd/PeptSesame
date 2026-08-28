"""
PeptSesame — Layer2: Multi-Evidence Scoring System
===================================================

Core innovation of the PeptSesame pipeline. Replaces experimental evidence
(Ribo-seq / mass spectrometry) with five computational evidence channels
for assessing the coding potential and biological relevance of small ORFs.

Evidence channels (all normalized to 0–1):
    A. Sequence features (k-mer composition, GC bias, hexamer scores)
    B. Cross-species conservation (BLASTP against related species peptides)
    C. RNA-seq expression evidence (if BAM/alignment files available)
    D. Structural features (signal peptide, transmembrane, cysteine-richness)
    E. Known motif / domain detection

Final score is a weighted average. Weights are configurable via YAML
or passed programmatically.

SLURM support:
    - Partition input into chunks, run EvidenceScorer per chunk
    - Merge results with merge_scores()

Typical usage::

    from pipeline.layer2_scoring.scoring_core import EvidenceScorer

    scorer = EvidenceScorer()
    result = scorer.evaluate(
        seq_id="sORF_0001",
        sequence="MASKLCYFFLFLFLVLLSLPSSHCDDDDDDDDDDDDDDDDDDE",
        orf_record={"start": 100, "end": 300, "frame": 1, "type": "intergenic"},
    )
    print(result["aggregated_score"])
"""

from __future__ import annotations

import gzip
import json
import itertools
import logging
import math
import os
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Standard k-mer lengths for compositional analysis
K_MER_SIZES = [1, 2, 3, 4]

# Default weight vector (production scoring, 4:3) — normalised to sum 1.0.
# The production ranking score uses two active core channels (sequence and
# structural features). Conservation, expression, and motif evidence are
# computed as independent orthogonal downstream layers (sub-scores retained
# in the output table) and are not included in the production ranking score.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "sequence_features": 4,
    "structural": 3,
}

# Amino acid groups for reduced-alphabet encoding
AA_HYDROPHOBIC = {"A", "V", "I", "L", "M", "F", "W", "Y", "P"}  # also includes Pro
AA_POLAR = {"S", "T", "C", "N", "Q"}
AA_POSITIVE = {"K", "R", "H"}
AA_NEGATIVE = {"D", "E"}
AA_SPECIAL = {"U", "O"}  # selenocysteine, pyrrolysine (rare)

# Valid amino acid one-letter codes
VALID_AA = frozenset(
    "ACDEFGHIKLMNPQRSTVWYXBZU*"
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EvidenceScores:
    """Container for all evidence sub-scores and the aggregate."""

    seq_id: str
    sequence_features: float = 0.0
    ml: float = 0.0                      # ML-based coding potential
    conservation: float = 0.0
    expression: float = 0.0
    structural: float = 0.0
    motif: float = 0.0
    aggregated_score: float = 0.0
    confidence: str = "low"  # "high" | "medium" | "low"

    def __post_init__(self) -> None:
        if not 0 <= self.sequence_features <= 1:
            raise ValueError("sequence_features must be in [0, 1]")
        if not 0 <= self.ml <= 1:
            raise ValueError("ml must be in [0, 1]")
        if not 0 <= self.conservation <= 1:
            raise ValueError("conservation must be in [0, 1]")
        if not 0 <= self.expression <= 1:
            raise ValueError("expression must be in [0, 1]")
        if not 0 <= self.structural <= 1:
            raise ValueError("structural must be in [0, 1]")
        if not 0 <= self.motif <= 1:
            raise ValueError("motif must be in [0, 1]")

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict (JSON-compatible)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "EvidenceScores":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ---------------------------------------------------------------------------
# Amino-acid helpers
# ---------------------------------------------------------------------------

def _validate_sequence(seq: str) -> str:
    """Strip whitespace, uppercase, and check for invalid characters."""
    seq = seq.strip().upper()
    invalid = set(seq) - VALID_AA
    if invalid:
        raise ValueError(
            f"Invalid amino acid characters in sequence: {''.join(sorted(invalid))}"
        )
    return seq


def _kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    """Compute normalized k-mer frequencies for a protein sequence.

    Returns a dict mapping each observed k-mer to its fraction of total k-mers.
    """
    if len(seq) < k:
        return {}
    total = len(seq) - k + 1
    freq: Dict[str, float] = {}
    for i in range(total):
        kmer = seq[i : i + k]
        freq[kmer] = freq.get(kmer, 0.0) + 1.0
    for kmer in freq:
        freq[kmer] /= total
    return freq


# ---------------------------------------------------------------------------
# Sequence feature computer
# ---------------------------------------------------------------------------

class SequenceFeatureComputer:
    """Compute compositional and statistical features from a protein sequence.

    These features feed into the coding-potential evidence score (Evidence A).
    """

    def __init__(self) -> None:
        self._hexamer_logprobs: Optional[Dict[str, float]] = None

    # ---- k-mer composition ----

    @staticmethod
    def compute_gc_content(seq: str) -> float:
        """Compute GC content of the nucleotide sequence (0-1).

        If a protein sequence is passed (amino acids), returns 0.0
        and logs a warning.
        """
        seq = seq.upper()
        # Detect if this looks like protein (any non-ACGT character)
        if set(seq) - {"A", "C", "G", "T", "N", "U"}:
            logger.warning(
                "compute_gc_content: sequence appears to be protein, "
                "returning 0.0"
            )
            return 0.0
        gc = seq.count("G") + seq.count("C")
        total = len(seq) - seq.count("N")
        return gc / total if total > 0 else 0.0

    @staticmethod
    def compute_amino_acid_composition(seq: str) -> Dict[str, float]:
        """Compute fractional composition of each amino acid (0-1)."""
        seq = _validate_sequence(seq)
        n = len(seq)
        if n == 0:
            return {}
        comp: Dict[str, float] = {}
        for aa in seq:
            comp[aa] = comp.get(aa, 0.0) + 1.0
        for aa in comp:
            comp[aa] /= n
        return comp

    @staticmethod
    def compute_reduced_alphabet_freqs(seq: str) -> Dict[str, float]:
        """Fraction of residues in each physico-chemical group."""
        seq = _validate_sequence(seq)
        n = len(seq)
        if n == 0:
            return {}
        groups = {
            "hydrophobic": AA_HYDROPHOBIC,
            "polar": AA_POLAR,
            "positive": AA_POSITIVE,
            "negative": AA_NEGATIVE,
        }
        counts: Dict[str, float] = {g: 0.0 for g in groups}
        for aa in seq:
            for gname, aa_set in groups.items():
                if aa in aa_set:
                    counts[gname] += 1.0
                    break
        return {g: c / n for g, c in counts.items()}

    @staticmethod
    def compute_kmer_freqs(seq: str, k: int = 2) -> Dict[str, float]:
        """Compute normalized di-/tri-/tetra-peptide frequencies."""
        return _kmer_frequencies(_validate_sequence(seq), k)

    @staticmethod
    def compute_kmer_diversity(seq: str, k: int = 2) -> float:
        """Shannon entropy of k-mer distribution — high = diverse.

        Ranges from 0 (single k-mer repeated) to log2(alphabet^k).
        """
        freqs = _kmer_frequencies(_validate_sequence(seq), k)
        if not freqs:
            return 0.0
        entropy = -sum(p * math.log2(p) for p in freqs.values())
        return entropy

    # ---- Hexamer score (coding potential) ----

    def load_hexamer_model(self, path: str) -> None:
        """Load pre-computed hexamer log-probabilities from JSON.

        Expected format: {"AAAAAA": -2.3, "AAAAAC": -1.8, ...}
        Can be trained from known coding / non-coding ORF hexamer counts.
        """
        with open(path) as f:
            self._hexamer_logprobs = json.load(f)
        logger.info("Loaded hexamer model with %d entries", len(self._hexamer_logprobs))

    def compute_hexamer_score(self, seq: str) -> float:
        """Compute average hexamer log-likelihood as coding-potential signal.

        Returns a value in 0-1 (sigmoid-normalised). Requires pre-loaded
        hexamer model. Returns 0.5 (neutral) if no model loaded.
        """
        if self._hexamer_logprobs is None:
            logger.debug("No hexamer model loaded — returning neutral 0.5")
            return 0.5
        seq = _validate_sequence(seq)
        if len(seq) < 6:
            return 0.5
        scores = []
        for i in range(len(seq) - 5):
            hexamer = seq[i : i + 6]
            score = self._hexamer_logprobs.get(hexamer, -3.0)  # default low
            scores.append(score)
        if not scores:
            return 0.5
        mean_lp = np.mean(scores)
        # Sigmoid with centre at 0 and scale 2
        return 1.0 / (1.0 + math.exp(-mean_lp / 2.0))

    # ---- Aggregate sequence-feature score ----

    def compute_score(self, seq: str) -> Dict[str, float]:
        """Compute a comprehensive sequence-feature score in [0, 1].

        Returns a dict with sub-scores and the final combined score.
        """
        seq = _validate_sequence(seq)
        length = len(seq)

        # 1. k-mer diversity (dipeptide diversity correlates with coding potential)
        dipep_div = self.compute_kmer_diversity(seq, k=2)
        max_dipep_entropy = math.log2(20 * 20)  # ~7.64
        norm_dipep = min(dipep_div / max_dipep_entropy, 1.0)

        # 2. Reduced-alphabet composition → measure of aa bias
        red = self.compute_reduced_alphabet_freqs(seq)
        # Hydrophobic fraction: typical coding seqs have moderate hydrophobicity
        hydro_score = 1.0 - abs(red.get("hydrophobic", 0) - 0.35) / 0.65

        # 3. Hexamer coding score (if model loaded)
        hex_score = self.compute_hexamer_score(seq)

        # 4. Length bonus: very short ORFs (<15 aa) are less likely to be coding
        length_bonus = min(length / 30.0, 1.0) if length > 0 else 0.0

        # 5. Amino acid variety (more distinct AAs → more coding-like)
        variety = len(set(seq)) / 20.0 if length > 0 else 0.0

        # Weighted combination
        combined = (
            0.25 * norm_dipep
            + 0.20 * hydro_score
            + 0.25 * hex_score
            + 0.15 * length_bonus
            + 0.15 * variety
        )

        return {
            "dipeptide_diversity": norm_dipep,
            "hydrophobic_balance": hydro_score,
            "hexamer_score": hex_score,
            "length_bonus": length_bonus,
            "aa_variety": variety,
            "combined": round(min(combined, 1.0), 4),
        }


# ---------------------------------------------------------------------------
# Conservation computer
# ---------------------------------------------------------------------------

class ConservationComputer:
    """Compute cross-species conservation evidence (Evidence B).

    Operates in two modes:
    - Live mode: runs BLASTP against a curated peptide database
    - Placeholder mode: returns 0.5 when BLAST DB not available
    """

    def __init__(self, blast_db: Optional[str] = None) -> None:
        self.blast_db = blast_db
        self._evalue_threshold = 1e-3
        self._identity_threshold = 30.0  # percent

    @property
    def available(self) -> bool:
        """Whether BLAST resources are available."""
        if self.blast_db is None:
            return False
        return os.path.isfile(self.blast_db) or os.path.isfile(self.blast_db + ".phr")

    def compute_score(
        self,
        seq_id: str,
        sequence: str,
        blast_results: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, float]:
        """Compute conservation score in [0, 1].

        Parameters
        ----------
        seq_id : str
            Identifier for logging / debugging.
        sequence : str
            Amino acid sequence.
        blast_results : list of dict, optional
            Pre-computed BLASTP hits. Each dict should have:
            - 'pident' : float — percent identity
            - 'evalue' : float — e-value
            - 'bitscore' : float — bit score
            - 'qcovs' : float — query coverage (percent)

        Returns
        -------
        dict with 'combined' key (the normalised score).
        """
        if blast_results is not None and len(blast_results) > 0:
            return self._score_from_blast(blast_results)
        if self.available:
            # Live BLAST (place holder — real call would subprocess)
            logger.info("BLAST DB available but not yet invoked for %s", seq_id)
            return {"combined": 0.5, "note": "blast_not_run"}
        return {"combined": 0.5, "note": "no_blast_db"}

    def _score_from_blast(
        self, hits: List[Dict[str, Any]]
    ) -> Dict[str, float]:
        """Score from a list of BLAST hits."""
        if not hits:
            return {"combined": 0.0, "n_hits": 0}

        # Filter by threshold
        good = [
            h
            for h in hits
            if h.get("evalue", 1) <= self._evalue_threshold
            and h.get("pident", 0) >= self._identity_threshold
        ]
        if not good:
            return {"combined": 0.0, "n_hits": 0}

        best = max(good, key=lambda h: h.get("bitscore", 0))
        pident = best.get("pident", 0) / 100.0
        qcov = best.get("qcovs", 100) / 100.0

        # Score: product of identity × coverage, sigmoid-shaped
        raw = pident * qcov
        score = min(raw * 1.5, 1.0)  # boost moderate hits slightly

        return {
            "combined": round(min(score, 1.0), 4),
            "n_hits": len(good),
            "best_pident": best.get("pident", 0),
            "best_evalue": best.get("evalue", 1),
            "best_bitscore": best.get("bitscore", 0),
        }


# ---------------------------------------------------------------------------
# Structural feature computer
# ---------------------------------------------------------------------------

class StructuralFeatureComputer:
    """Compute structural evidence score (Evidence D).

    Includes:
    - Signal peptide prediction (rule-based or via external tool)
    - Transmembrane domain detection
    - Cysteine-rich region analysis (CRP)
    - Disulfide bond pattern
    """

    # Eukaryotic signal peptide heuristics
    # SignalP-like: N-region (basic), H-region (hydrophobic), C-region (neutral)
    SIGNALP_CUTOFF = 0.5  # placeholder for real SignalP integration

    def __init__(self) -> None:
        self.signalp_path: Optional[str] = None

    # ---- Signal peptide ----

    @staticmethod
    def predict_signal_peptide_rule_based(seq: str) -> Tuple[float, Optional[str]]:
        """Rule-based signal peptide prediction.

        Returns (score, cleavage_site) where score is in [0, 1].
        Simple heuristic: checks for N-terminal hydrophobic stretch
        followed by a cleavage motif (A-X-A or similar).
        """
        seq = _validate_sequence(seq)
        if len(seq) < 20:
            return (0.0, None)

        # Look for hydrophobic core in first 20 residues
        first20 = seq[:20]
        hydrophobic_run = max(
            (len(list(g)) for c, g in itertools.groupby(first20, key=lambda aa: aa in AA_HYDROPHOBIC) if c),
            default=0,
        )

        if hydrophobic_run < 6:
            return (0.0, None)

        # Check for potential cleavage site: residues 15-25 typically
        # Consensus: small neutral residue at -1 and -3 (A, G, S, C, T)
        # Patterns: [AVI].[AG]| or similar
        for pos in range(14, min(24, len(seq) - 3)):
            tripep = seq[pos : pos + 3]
            # Simple heuristic: A/G/S at pos+2, small at pos
            if tripep[0] in "ASGT" and tripep[2] in "ASGC":
                score = min(0.5 + hydrophobic_run * 0.03, 1.0)
                return (round(score, 4), f"{pos + 2}-{pos + 3}")
        return (0.0, None)

    def predict_signal_peptide_signalp(self, seq: str) -> Tuple[float, Optional[str]]:
        """Integration point for SignalP (external tool).

        Placeholder — requires SignalP binary in PATH or configured path.
        Returns (0.5, None) if SignalP not found.
        """
        if self.signalp_path is None or not os.path.isfile(self.signalp_path):
            return (0.5, None)
        # Real implementation would subprocess check SignalP
        logger.debug("SignalP binary found at %s — not yet invoked", self.signalp_path)
        return (0.5, None)

    # ---- Transmembrane ----

    @staticmethod
    def predict_transmembrane_helices(seq: str) -> Tuple[float, int]:
        """Simple TM helix predictor based on hydrophobic stretches.

        Returns (score, n_tm_helices). Score in [0, 1].
        A single TM helix is normal — zero or 2+ reduces score slightly.
        """
        seq = _validate_sequence(seq)
        if len(seq) < 21:
            return (0.0, 0)

        # Sliding window of 21 aa: hydrophobic content
        n_tm = 0
        i = 0
        while i <= len(seq) - 21:
            window = seq[i : i + 21]
            hydro_count = sum(1 for aa in window if aa in AA_HYDROPHOBIC)
            if hydro_count >= 14:  # ≥ 2/3 hydrophobic
                n_tm += 1
                i += 21
            else:
                i += 1

        # Score logic:
        # 0 TM → low (signal peptide might be mis-called)
        # 1 TM → good (typical SSP)
        # 2+ → possible membrane protein, less likely as signalling peptide
        if n_tm == 0:
            score = 0.4
        elif n_tm == 1:
            score = 0.9
        elif n_tm == 2:
            score = 0.5
        else:
            score = 0.2

        return (round(score, 4), n_tm)

    # ---- Cysteine richness ----

    @staticmethod
    def compute_cysteine_rich_score(seq: str) -> float:
        """Score for cysteine-rich profile (CRP).

        CRPs have a distinct Cys distribution (typically 6 or 8 Cys).
        Returns score in [0, 1].
        """
        seq = _validate_sequence(seq)
        if len(seq) < 20:
            return 0.0
        n_cys = seq.count("C")
        cys_freq = n_cys / len(seq)

        # Typical CRP: 4–12% Cys
        ideal_low, ideal_high = 0.04, 0.12
        if ideal_low <= cys_freq <= ideal_high:
            return min(0.5 + n_cys * 0.05, 1.0)
        if cys_freq > ideal_high:
            return max(0.3, 1.0 - (cys_freq - ideal_high) * 2)
        return max(0.0, cys_freq / ideal_low * 0.5)

    @staticmethod
    def predict_disulfide_pattern(seq: str) -> int:
        """Count predicted disulfide bonds (pairs of Cys).

        Simple: assumes even-numbered Cys can form S-S bonds.
        Returns number of potential bonds.
        """
        cys_positions = [i for i, aa in enumerate(seq) if aa == "C"]
        return len(cys_positions) // 2

    # ---- Aggregate ----

    def compute_score(self, seq: str) -> Dict[str, float]:
        """Aggregate structural feature score."""
        sp_score, sp_site = self.predict_signal_peptide_rule_based(seq)
        tm_score, n_tm = self.predict_transmembrane_helices(seq)
        crp_score = self.compute_cysteine_rich_score(seq)

        # Weight: signal peptide highest, then TM, then CRP
        combined = 0.50 * sp_score + 0.30 * tm_score + 0.20 * crp_score
        combined = round(min(combined, 1.0), 4)

        return {
            "signal_peptide_score": sp_score,
            "signal_peptide_site": sp_site or "",
            "tm_score": tm_score,
            "n_tm_helices": n_tm,
            "crp_score": crp_score,
            "n_disulfide": self.predict_disulfide_pattern(seq),
            "combined": combined,
        }


# ---------------------------------------------------------------------------
# Expression computer
# ---------------------------------------------------------------------------

class ExpressionComputer:
    """Compute RNA-seq expression evidence (Evidence C).

    Operates in three modes:
    1. Direct: read from pre-computed expression table (tsv/csv)
    2. BAM-based: count reads overlapping sORF regions (requires pysam)
    3. Placeholder: returns neutral 0.5 when no data
    """

    def __init__(self) -> None:
        self._expression_data: Dict[str, float] = {}  # seq_id → TPM/RPKM

    def load_expression_table(self, path: str) -> None:
        """Load a pre-computed expression table.

        Expected format (TSV)::
            seq_id\texpression_value
            sORF_0001\t15.3
        """
        expr: Dict[str, float] = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("seq_id"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    expr[parts[0]] = float(parts[1])
        self._expression_data.update(expr)
        logger.info("Loaded expression data for %d entries", len(expr))

    def compute_score(
        self,
        seq_id: str,
        orf_record: Optional[Dict[str, Any]] = None,
        bam_path: Optional[str] = None,
    ) -> Dict[str, float]:
        """Compute expression evidence score.

        Parameters
        ----------
        seq_id : str
            sORF identifier.
        orf_record : dict, optional
            Genomic coordinates for BAM-based counting.
            Keys: 'chrom', 'start', 'end', 'strand'.
        bam_path : str, optional
            Path to indexed BAM file.

        Returns
        -------
        dict with 'combined' score and expression value.
        """
        # Mode 1: pre-computed table
        if seq_id in self._expression_data:
            rpkm = self._expression_data[seq_id]
            return self._score_from_rpkm(rpkm)

        # Mode 2: BAM-based (placeholder)
        if bam_path and orf_record:
            logger.debug("BAM-based expression not yet implemented for %s", seq_id)
            return {"combined": 0.5, "rpkm": 0.0, "note": "bam_not_processed"}

        # No data
        return {"combined": 0.5, "rpkm": 0.0, "note": "no_expression_data"}

    @staticmethod
    def _score_from_rpkm(rpkm: float) -> Dict[str, float]:
        """Convert RPKM/TPM value to a [0, 1] score.

        Thresholds (typical for plant sORFs):
        - RPKM < 0.5: essentially no expression
        - RPKM 0.5–2: weak
        - RPKM 2–10: moderate
        - RPKM > 10: strong
        """
        if rpkm <= 0:
            return {"combined": 0.0, "rpkm": rpkm}

        # Logistic function centered at 1.0, slope 0.5
        score = 1.0 / (1.0 + math.exp(-0.5 * (rpkm - 1.0)))
        # But cap min at 0.1 for any detectable expression
        if rpkm > 0.01:
            score = max(score, 0.1)

        return {"combined": round(min(score, 1.0), 4), "rpkm": round(rpkm, 4)}


# ---------------------------------------------------------------------------
# Motif / domain computer
# ---------------------------------------------------------------------------

class MotifComputer:
    """Compute known motif / domain evidence (Evidence E).

    Matches sequences against curated motif profiles for:
    - Secreted signalling peptide families (CLE, RALF, CEP, etc.)
    - Known functional domains (e.g., cysteine-knot, defensin-like)
    - General small-peptide motifs
    """

    def __init__(self, motif_db: Optional[Dict[str, str]] = None) -> None:
        self.motif_db = motif_db or {}  # name → regex pattern
        self._compile_motifs()

    def _compile_motifs(self) -> None:
        self._compiled: List[Tuple[str, re.Pattern]] = []
        for name, pattern in self.motif_db.items():
            try:
                self._compiled.append((name, re.compile(pattern)))
            except re.error as e:
                logger.warning("Invalid regex for motif '%s': %s", name, e)

    def load_motif_file(self, path: str) -> None:
        """Load motifs from a JSON file: {"motif_name": "regexp", ...}."""
        with open(path) as f:
            motifs = json.load(f)
        self.motif_db.update(motifs)
        self._compile_motifs()

    def compute_score(self, seq: str) -> Dict[str, float]:
        """Match sequence against known motif profiles.

        Returns combined score and list of matched motifs.
        """
        seq = _validate_sequence(seq)
        matches: List[str] = []
        for name, pattern in self._compiled:
            if pattern.search(seq):
                matches.append(name)

        if matches:
            # More matches → higher confidence, but saturate at 3+
            motif_count = min(len(matches), 5)
            score = 0.3 + 0.15 * motif_count
        else:
            score = 0.0

        return {
            "combined": round(min(score, 1.0), 4),
            "matched_motifs": matches,
        }


# ===========================================================================
# Evidence F: ML-based coding potential scorer
# ===========================================================================

class MlScorerComputer:
    """ML-based coding potential scorer — wraps LightGBMPredictor.

    Provides the same compute_score(seq_id, sequence) interface as other
    evidence computers, returning a dict with key "combined".

    Falls back gracefully to rule-based scoring when ML model is unavailable.
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        threshold_high: float = 0.7,
        threshold_medium: float = 0.4,
    ):
        self.model_path = model_path
        self._predictor = None

    def _get_predictor(self):
        if self._predictor is None:
            from .ml_scorer import LightGBMPredictor
            if self.model_path and os.path.exists(self.model_path):
                self._predictor = LightGBMPredictor(model_path=self.model_path)
            else:
                self._predictor = LightGBMPredictor(model=None)
        return self._predictor

    def compute_score(self, seq_id: str, sequence: str) -> Dict:
        """Compute ML-based coding potential score.

        Returns dict with 'combined' key (0-1) matching other evidence formats.
        """
        predictor = self._get_predictor()
        result = predictor.predict(sequence)

        return {
            "combined": result["ml_score"],
            "ml_confidence": result["confidence"],
            "ml_probability": result["probability"],
            "model_ready": predictor.is_ready(),
        }

    def compute_score_batch(self, records: List[Dict]) -> List[Dict]:
        """Batch scoring for efficiency."""
        predictor = self._get_predictor()
        seqs = [r["sequence"] for r in records]
        results = predictor.predict_batch(seqs)

        out = []
        for r, res in zip(records, results):
            out.append({
                "combined": res["ml_score"],
                "ml_confidence": res["confidence"],
                "ml_probability": res["probability"],
            })
        return out


# ---------------------------------------------------------------------------
# Main scorer
# ---------------------------------------------------------------------------

class EvidenceScorer:
    """Multi-evidence scoring system for small ORF coding potential.

    Combines the active core channels into a single aggregated ranking score.
    Supports configurable weights, per-channel sub-scores, and
    SLURM-aware batch processing.

    Parameters
    ----------
    weights : dict, optional
        Evidence weights. Keys: 'sequence_features', 'conservation',
        'expression', 'structural', 'motif'. Defaults to DEFAULT_WEIGHTS.
    motif_db : dict, optional
        Motif regex patterns.
    blast_db : str, optional
        Path to BLAST database for conservation search.
    expression_table : str, optional
        Path to pre-computed expression data (TSV).
    signalp_path : str, optional
        Path to SignalP binary.
    """

    def __init__(
        self,
        weights: Optional[Dict[str, float]] = None,
        motif_db: Optional[Dict[str, str]] = None,
        blast_db: Optional[str] = None,
        expression_table: Optional[str] = None,
        signalp_path: Optional[str] = None,
        ml_model_path: Optional[str] = None,
    ) -> None:
        self.weights = weights or dict(DEFAULT_WEIGHTS)

        # Normalise weights
        total = sum(self.weights.values())
        if abs(total - 1.0) > 0.01:
            logger.warning(
                "Weights sum to %.3f (not 1.0) — normalising", total
            )
            self.weights = {k: v / total for k, v in self.weights.items()}

        # Evidence computers (lazy init)
        self._seq_feat: Optional[SequenceFeatureComputer] = None
        self._conservation: Optional[ConservationComputer] = None
        self._structural: Optional[StructuralFeatureComputer] = None
        self._expression: Optional[ExpressionComputer] = None
        self._motif: Optional[MotifComputer] = None
        self._ml_scorer: Optional["MlScorerComputer"] = None

        # Pass-through config
        self._blast_db = blast_db
        self._expression_table = expression_table
        self._signalp_path = signalp_path
        self._motif_db = motif_db
        self._ml_model_path = ml_model_path

    # ---- Lazy accessors ----

    @property
    def seq_feat(self) -> SequenceFeatureComputer:
        if self._seq_feat is None:
            self._seq_feat = SequenceFeatureComputer()
        return self._seq_feat

    @property
    def conservation(self) -> ConservationComputer:
        if self._conservation is None:
            self._conservation = ConservationComputer(blast_db=self._blast_db)
        return self._conservation

    @property
    def structural(self) -> StructuralFeatureComputer:
        if self._structural is None:
            self._structural = StructuralFeatureComputer()
            self._structural.signalp_path = self._signalp_path
        return self._structural

    @property
    def expression(self) -> ExpressionComputer:
        if self._expression is None:
            comp = ExpressionComputer()
            if self._expression_table:
                comp.load_expression_table(self._expression_table)
            self._expression = comp
        return self._expression

    @property
    def motif(self) -> MotifComputer:
        if self._motif is None:
            self._motif = MotifComputer(motif_db=self._motif_db)
        return self._motif

    @property
    def ml_scorer(self) -> "MlScorerComputer":
        if self._ml_scorer is None:
            self._ml_scorer = MlScorerComputer(model_path=self._ml_model_path)
        return self._ml_scorer

    # ---- Public API ----

    def evaluate(
        self,
        seq_id: str,
        sequence: str,
        orf_record: Optional[Dict[str, Any]] = None,
        blast_results: Optional[List[Dict[str, Any]]] = None,
        bam_path: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run all evidence channels and return the aggregated result.

        Parameters
        ----------
        seq_id : str
            Unique identifier for this sORF.
        sequence : str
            Amino acid sequence (single letter code).
        orf_record : dict, optional
            Genomic context: {'chrom', 'start', 'end', 'strand', 'type'}.
        blast_results : list of dict, optional
            Pre-computed BLASTP hits.
        bam_path : str, optional
            Path to indexed BAM file for expression counting.

        Returns
        -------
        dict
            Full evaluation result with per-channel scores and aggregate.
        """
        seq = _validate_sequence(sequence)

        # Evidence A: sequence features
        feat = self.seq_feat.compute_score(seq)

        # Evidence B: conservation
        cons = self.conservation.compute_score(seq_id, seq, blast_results)

        # Evidence C: expression
        exp = self.expression.compute_score(seq_id, orf_record, bam_path)

        # Evidence D: structural
        struct = self.structural.compute_score(seq)

        # Evidence E: motif
        mot = self.motif.compute_score(seq)

        # Evidence F: ML-based coding potential
        ml = self.ml_scorer.compute_score(seq_id, seq)

        # Aggregate
        scores = EvidenceScores(
            seq_id=seq_id,
            sequence_features=feat["combined"],
            ml=ml["combined"],
            conservation=cons["combined"],
            expression=exp["combined"],
            structural=struct["combined"],
            motif=mot["combined"],
        )

        aggregated = self.aggregate_score(scores)
        scores_dict = scores.to_dict()
        scores_dict["aggregated_score"] = aggregated

        # Confidence tiers
        if aggregated >= 0.7:
            scores_dict["confidence"] = "high"
        elif aggregated >= 0.4:
            scores_dict["confidence"] = "medium"
        else:
            scores_dict["confidence"] = "low"

        # Attach sub-scores for transparency
        scores_dict["sub_scores"] = {
            "sequence_features": feat,
            "conservation": cons,
            "expression": exp,
            "structural": struct,
            "motif": mot,
            "ml": ml,
        }

        return scores_dict

    def aggregate_score(
        self, scores: EvidenceScores
    ) -> float:
        """Weighted sum of the active core channel scores.

        Weights are taken from ``self.weights``.
        """
        raw = (
            self.weights["sequence_features"] * scores.sequence_features
            + self.weights["ml"] * scores.ml
            + self.weights["conservation"] * scores.conservation
            + self.weights["expression"] * scores.expression
            + self.weights["structural"] * scores.structural
            + self.weights["motif"] * scores.motif
        )
        return round(min(raw, 1.0), 4)

    # ---- Batch processing ----

    def evaluate_batch(
        self,
        records: List[Dict[str, Any]],
        batch_size: int = 10000,
        n_jobs: int = 1,
    ) -> List[Dict[str, Any]]:
        """Evaluate multiple sORFs, optionally in parallel.

        Parameters
        ----------
        records : list of dict
            Each dict must have 'seq_id' and 'sequence'.
            May also have 'orf_record', 'blast_results', 'bam_path'.
        batch_size : int
            For future SLURM chunking — size of each chunk.
        n_jobs : int
            Number of parallel workers (1 = serial).

        Returns
        -------
        list of dict
            Evaluation results for each record.
        """
        results: List[Dict[str, Any]] = []

        if n_jobs <= 1:
            for rec in records:
                res = self.evaluate(
                    seq_id=rec["seq_id"],
                    sequence=rec["sequence"],
                    orf_record=rec.get("orf_record"),
                    blast_results=rec.get("blast_results"),
                    bam_path=rec.get("bam_path"),
                )
                results.append(res)
        else:
            # Placeholder for multiprocessing / SLURM
            logger.warning(
                "Parallel evaluation (n_jobs=%d) not yet implemented — "
                "falling back to serial",
                n_jobs,
            )
            for rec in records:
                res = self.evaluate(
                    seq_id=rec["seq_id"],
                    sequence=rec["sequence"],
                    orf_record=rec.get("orf_record"),
                    blast_results=rec.get("blast_results"),
                    bam_path=rec.get("bam_path"),
                )
                results.append(res)

        return results

    # ---- SLURM support ----

    @staticmethod
    def chunk_for_slurm(
        records: List[Dict[str, Any]],
        n_chunks: int,
        output_dir: str,
        prefix: str = "scoring_chunk",
    ) -> List[str]:
        """Split records into chunks for SLURM array jobs.

        Parameters
        ----------
        records : list of dict
            Input sORF records.
        n_chunks : int
            Number of SLURM array tasks.
        output_dir : str
            Directory to write chunk files.
        prefix : str
            Prefix for chunk filenames.

        Returns
        -------
        list of str
            Paths to written chunk JSON files.
        """
        os.makedirs(output_dir, exist_ok=True)
        chunk_paths: List[str] = []

        for i in range(n_chunks):
            chunk = records[i::n_chunks]
            path = os.path.join(output_dir, f"{prefix}_{i:04d}.json")
            with open(path, "w") as f:
                json.dump(chunk, f, indent=2)
            chunk_paths.append(path)
            logger.info("Wrote %d records to %s", len(chunk), path)

        return chunk_paths

    @staticmethod
    def merge_scores(
        chunk_dir: str,
        output_path: str,
        pattern: str = "scoring_chunk_*.json.results",
    ) -> str:
        """Merge per-chunk scoring results into a single JSON file.

        For use as the final step in a SLURM array pipeline.
        """
        import glob

        merged: List[Dict[str, Any]] = []
        for fpath in sorted(glob.glob(os.path.join(chunk_dir, pattern))):
            with open(fpath) as f:
                data = json.load(f)
                if isinstance(data, list):
                    merged.extend(data)
                else:
                    merged.append(data)

        with open(output_path, "w") as f:
            json.dump(merged, f, indent=2)

        logger.info(
            "Merged scores from %s → %s (%d records)",
            chunk_dir, output_path, len(merged),
        )
        return output_path


# ---------------------------------------------------------------------------
# Convenience CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Minimal CLI for testing / debugging."""
    import argparse

    parser = argparse.ArgumentParser(
        description="PeptSesame Layer2 — Multi-evidence sORF scoring"
    )
    parser.add_argument("--sequence", "-s", required=True, help="Amino acid sequence")
    parser.add_argument("--seq-id", "-i", default="test", help="Sequence identifier")
    parser.add_argument(
        "--weights", nargs=5, type=float,
        metavar=("SEQ", "CON", "EXP", "STR", "MOT"),
        default=[0.30, 0.20, 0.15, 0.20, 0.15],
        help="Evidence weights (5 floats)",
    )
    args = parser.parse_args()

    scorer = EvidenceScorer(
        weights={
            "sequence_features": args.weights[0],
            "conservation": args.weights[1],
            "expression": args.weights[2],
            "structural": args.weights[3],
            "motif": args.weights[4],
        }
    )

    result = scorer.evaluate(args.seq_id, args.sequence)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
