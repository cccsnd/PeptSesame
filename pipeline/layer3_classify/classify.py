"""
PeptSesame — Layer3: Small Peptide Classification
==================================================

Classifies candidate sORF-encoded peptides into functional categories:

    1. Classic SSP families   — CLE, RALF, CEP, PSK, PSY1, IDA, EPFL, RGF
    2. miPEP                  — miRNA-encoded peptides (pri-miRNA sORFs)
    3. uORF                   — upstream ORFs (5' UTR located)
    4. lncORF                 — lncRNA-encoded peptides
    5. AMP                    — antimicrobial peptides (database comparison)

Classification is multi-pass: a single sORF may match multiple categories
(e.g., an SSP-like peptide in a lncRNA). Each assignment carries a
confidence score (0-1).

SLURM support:
    - Chunk input GFF/FASTA for parallel classification
    - Merge per-chunk reports with merge_classifications()

Typical usage::

    from pipeline.layer3_classify.classify import SmallPeptideClassifier

    classifier = SmallPeptideClassifier()
    result = classifier.classify_sorf(
        seq_id="sORF_0001",
        sequence="MASKLCYFFLFLFLVLLSLPSSH...",
        gff_record={"chrom": "Chr1", "start": 1000, "end": 1300,
                     "strand": "+", "feature_type": "CDS"},
    )
    print(result["categories"])
"""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import logging
import math
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from pipeline.layer3_classify.motif_profiles import (
    SSP_MOTIFS,
    SSP_ALIASES,
    SSP_FAMILY_INFO,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Classification categories
CATEGORY_SSP = "SSP"          # Small Signalling Peptide (CLE, RALF, etc.)
CATEGORY_MIPEP = "miPEP"      # miRNA-encoded peptide
CATEGORY_UORF = "uORF"        # Upstream ORF (5' UTR)
CATEGORY_LNCORF = "lncORF"    # lncRNA-encoded ORF
CATEGORY_AMP = "AMP"          # Antimicrobial peptide

ALL_CATEGORIES = [
    CATEGORY_SSP,
    CATEGORY_MIPEP,
    CATEGORY_UORF,
    CATEGORY_LNCORF,
    CATEGORY_AMP,
]

# Minimum ORF length (aa) for each category
MIN_LENGTH: Dict[str, int] = {
    CATEGORY_SSP: 20,
    CATEGORY_MIPEP: 8,
    CATEGORY_UORF: 5,
    CATEGORY_LNCORF: 5,
    CATEGORY_AMP: 10,
}

# Maximum ORF length (aa) for each category
MAX_LENGTH: Dict[str, int] = {
    CATEGORY_SSP: 300,
    CATEGORY_MIPEP: 100,
    CATEGORY_UORF: 100,
    CATEGORY_LNCORF: 150,
    CATEGORY_AMP: 150,
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ClassificationResult:
    """Result of classifying a single sORF."""

    seq_id: str
    sequence: str
    length_aa: int
    categories: Dict[str, float] = field(default_factory=dict)
    # category_name → confidence (0-1)
    ssp_family: Optional[str] = None
    # If SSP, which family (CLE, RALF, etc.)
    location_context: str = "intergenic"
    # intergenic / utr5 / cds / lncrna / intron / mirna
    notes: List[str] = field(default_factory=list)
    amp_match: Optional[Dict[str, Any]] = None
    # Details of AMP database match

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ClassificationResult":
        """Deserialize from dict."""
        return cls(**{k: v for k, v in d.items()
                      if k in cls.__dataclass_fields__})


# ---------------------------------------------------------------------------
# SSP Family classifier
# ---------------------------------------------------------------------------

class SSPFamilyClassifier:
    """Classify sORFs into known signalling peptide families by motif matching.

    Uses curated regex patterns from motif_profiles.py.
    Supports both full-length and active-domain matching.
    """

    def __init__(self, custom_motifs: Optional[Dict[str, str]] = None) -> None:
        motifs = dict(SSP_MOTIFS)
        if custom_motifs:
            motifs.update(custom_motifs)

        self._patterns: Dict[str, re.Pattern] = {}
        for family, pattern in motifs.items():
            try:
                self._patterns[family] = re.compile(pattern, re.IGNORECASE)
            except re.error as e:
                logger.warning("Invalid regex for SSP family '%s': %s", family, e)

    def classify(self, seq: str) -> Dict[str, Any]:
        """Match a sequence against SSP family motifs.

        Returns dict with:
        - 'ssp_families': list of (family_name, pattern_match_detail)
        - 'confidence': best match confidence (0-1)
        - 'best_family': highest-scoring family (or None)
        """
        seq = seq.strip().upper()
        matches: List[Tuple[str, float, str]] = []  # (family, confidence, detail)

        for family, pattern in self._patterns.items():
            m = pattern.search(seq)
            if m:
                matched_seq = m.group()
                # Confidence based on match length relative to typical active peptide
                match_len = len(matched_seq)
                # Longer matches → more confident
                confidence = min(0.4 + match_len * 0.02, 0.95)
                matches.append((family, round(confidence, 4), matched_seq))

        if not matches:
            return {
                "ssp_families": [],
                "confidence": 0.0,
                "best_family": None,
            }

        matches.sort(key=lambda x: x[1], reverse=True)

        return {
            "ssp_families": [
                {"family": fam, "confidence": conf, "matched_region": mseq}
                for fam, conf, mseq in matches
            ],
            "confidence": matches[0][1],
            "best_family": matches[0][0],
        }


# ---------------------------------------------------------------------------
# Location-based classifier
# ---------------------------------------------------------------------------

class LocationClassifier:
    """Classify sORFs by their genomic location context.

    Uses GFF annotations to determine whether an sORF falls within:
    - 5' UTR (→ uORF)
    - CDS (→ possible alternative ORF, excluded by Layer1)
    - lncRNA transcript (→ lncORF)
    - Intron
    - miRNA precursor (→ miPEP)
    - Intergenic region
    """

    def __init__(self, gff_path: Optional[str] = None) -> None:
        self._gff_index: Optional[Dict[str, List[Dict[str, Any]]]] = None
        if gff_path:
            self.load_gff(gff_path)

    def load_gff(self, path: str) -> None:
        """Load and index a GFF/GTF annotation file.

        Indexes features by chromosome for fast overlap queries.
        Supports both GFF3 and GTF formats (simple parser).
        """
        logger.info("Loading GFF annotations from %s", path)
        index: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                chrom, source, ftype, start, end, score, strand, phase, attr = parts
                try:
                    start_i = int(start)
                    end_i = int(end)
                except ValueError:
                    continue
                record = {
                    "chrom": chrom,
                    "source": source,
                    "type": ftype,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "attributes": attr,
                }
                index[chrom].append(record)

        self._gff_index = dict(index)
        logger.info(
            "Loaded %d features across %d chromosomes",
            sum(len(v) for v in index.values()),
            len(index),
        )

    def classify_location(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: str,
    ) -> Dict[str, Any]:
        """Determine the genomic context of an sORF interval.

        Returns dict with:
        - 'context': one of 'utr5', 'cds', 'lncrna', 'intron', 'mirna', 'intergenic'
        - 'overlapping_features': names of overlapping annotations
        - 'confidence': confidence in the assignment (0-1)
        """
        if self._gff_index is None:
            return {
                "context": "intergenic",
                "overlapping_features": [],
                "confidence": 0.5,
                "note": "no_gff_loaded",
            }

        features = self._gff_index.get(chrom, [])
        overlapping: List[Dict[str, Any]] = []
        for feat in features:
            if feat["end"] >= start and feat["start"] <= end:
                overlapping.append(feat)

        if not overlapping:
            return {
                "context": "intergenic",
                "overlapping_features": [],
                "confidence": 0.8,
            }

        # Priority order: miRNA > lncRNA > 5'UTR > CDS > intron
        context_order = ["miRNA", "lnc_RNA", "lncRNA", "five_prime_UTR",
                         "UTR5", "5utr", "CDS", "intron", "exon"]

        for ctx_type in context_order:
            for feat in overlapping:
                ftype = feat["type"].lower().replace(" ", "_")
                if ctx_type.lower() in ftype:
                    mapped = self._map_type(ftype)
                    return {
                        "context": mapped,
                        "overlapping_features": [
                            {"type": feat["type"],
                             "start": feat["start"],
                             "end": feat["end"],
                             "strand": feat["strand"]}
                        ],
                        "confidence": 0.9,
                    }

        # Check for any transcript feature
        for feat in overlapping:
            ftype = feat["type"].lower()
            if "transcript" in ftype or "mrna" in ftype:
                return {
                    "context": "intergenic",
                    "overlapping_features": [{"type": feat["type"],
                                               "start": feat["start"],
                                               "end": feat["end"],
                                               "strand": feat["strand"]}],
                    "confidence": 0.5,
                    "note": "overlaps_transcript_no_feature_type",
                }

        return {
            "context": "intergenic",
            "overlapping_features": [{"type": f["type"],
                                       "start": f["start"],
                                       "end": f["end"]} for f in overlapping[:3]],
            "confidence": 0.3,
        }

    @staticmethod
    def _map_type(ftype: str) -> str:
        """Map GFF feature type to classification category."""
        if "mirna" in ftype:
            return "mirna"
        if "lnc" in ftype:
            return "lncrna"
        if "utr5" in ftype or "five_prime" in ftype:
            return "utr5"
        if "cds" in ftype:
            return "cds"
        if "intron" in ftype:
            return "intron"
        return "intergenic"

    # ---- Direct context assignment ----

    @staticmethod
    def classify_by_context(
        gff_record: Dict[str, Any],
    ) -> str:
        """Classify based on a pre-parsed GFF record.

        Key in gff_record: 'feature_type' or 'type'.
        Returns one of: 'utr5', 'cds', 'lncrna', 'intron', 'mirna', 'intergenic'.
        """
        ftype = gff_record.get("feature_type", gff_record.get("type", "")).lower()
        if "mirna" in ftype or ftype == "miRNA":
            return "mirna"
        if "lnc" in ftype:
            return "lncrna"
        if "utr5" in ftype or ftype == "five_prime_UTR":
            return "utr5"
        if "cds" in ftype:
            return "cds"
        if "intron" in ftype:
            return "intron"
        return "intergenic"


# ---------------------------------------------------------------------------
# AMP classifier
# ---------------------------------------------------------------------------

class AMPClassifier:
    """Classify sORFs as antimicrobial peptides (AMPs).

    Uses two strategies:
    1. Sequence similarity to known AMP databases
    2. Physico-chemical properties (cationic, amphipathic, small)

    Known AMP databases are loaded as FASTA or JSON.
    """

    def __init__(self, amp_db_path: Optional[str] = None) -> None:
        self._amp_sequences: List[Tuple[str, str]] = []  # (id, sequence)
        if amp_db_path:
            self.load_amp_database(amp_db_path)

    def load_amp_database(self, path: str) -> None:
        """Load AMP sequences from FASTA file."""
        logger.info("Loading AMP database from %s", path)
        seq_id = ""
        seq_chars: List[str] = []
        count = 0
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq_id and seq_chars:
                        self._amp_sequences.append((seq_id, "".join(seq_chars)))
                        count += 1
                    seq_id = line[1:].split()[0]
                    seq_chars = []
                elif line:
                    seq_chars.append(line.upper())
        if seq_id and seq_chars:
            self._amp_sequences.append((seq_id, "".join(seq_chars)))
            count += 1
        logger.info("Loaded %d AMP sequences", count)

    # ---- Physico-chemical AMP prediction ----

    @staticmethod
    def predict_amp_properties(seq: str) -> Dict[str, float]:
        """Predict AMP potential from sequence properties.

        Returns dict with property scores and a combined 'score'.
        Uses the classic AMP criteria:
        - Net positive charge (+2 to +9)
        - High isoelectric point (pI > 8)
        - Short length (10-50 aa)
        - Hydrophobic content 30-50%
        - Amphipathic structure indicator
        """
        seq = seq.strip().upper()
        n = len(seq)
        if n == 0:
            return {"score": 0.0, "isoelectric_point": 0.0,
                    "net_charge": 0.0, "hydrophobic_ratio": 0.0}

        # Net charge at pH 7 (approximate)
        # K=+1, R=+1, H=+0.5, D=-1, E=-1
        net_charge = (
            seq.count("K") + seq.count("R") + 0.5 * seq.count("H")
            - seq.count("D") - seq.count("E")
        )

        # pI approximation (very rough)
        pos = seq.count("K") + seq.count("R")
        neg = seq.count("D") + seq.count("E")
        pi = 7.0 + 2.0 * (pos - neg) / max(n, 1)

        # Hydrophobic ratio
        hydro_count = sum(1 for aa in seq if aa in "AILMFWVG")
        hydro_ratio = hydro_count / n

        # Length score (AMPs typically 10-50 aa)
        len_score = 1.0 - abs(n - 30) / 50.0 if n >= 10 else 0.0
        len_score = max(0.0, min(1.0, len_score))

        # Charge score
        charge_score = 0.0
        if 2 <= net_charge <= 9:
            charge_score = 0.8 + 0.02 * (net_charge - 2)
        elif net_charge > 9:
            charge_score = 0.9
        elif net_charge > 0:
            charge_score = 0.3
        charge_score = min(charge_score, 1.0)

        # Hydrophobicity score (ideal: 30-50%)
        hydro_score = 1.0 - abs(hydro_ratio - 0.4) / 0.6
        hydro_score = max(0.0, min(1.0, hydro_score))

        # pI score (AMPs typically basic, pI > 8)
        pi_score = min((pi - 7.0) / 5.0, 1.0) if pi > 7.0 else 0.0
        pi_score = max(0.0, pi_score)

        combined = 0.30 * charge_score + 0.25 * hydro_score + \
                   0.25 * len_score + 0.20 * pi_score
        combined = round(min(combined, 1.0), 4)

        return {
            "score": combined,
            "net_charge": round(net_charge, 1),
            "isoelectric_point": round(pi, 1),
            "hydrophobic_ratio": round(hydro_ratio, 4),
            "length_score": round(len_score, 4),
            "charge_score": round(charge_score, 4),
            "hydro_score": round(hydro_score, 4),
            "pi_score": round(pi_score, 4),
        }

    def classify(self, seq: str) -> Dict[str, Any]:
        """Classify a sequence as AMP or not.

        Returns dict with:
        - 'is_amp': whether it matches AMP criteria
        - 'confidence': confidence score (0-1)
        - 'properties': physico-chemical property scores
        - 'best_db_match': best database hit, if any
        """
        seq = seq.strip().upper()
        n = len(seq)

        result: Dict[str, Any] = {
            "is_amp": False,
            "confidence": 0.0,
            "properties": {},
            "best_db_match": None,
        }

        # 1. Physico-chemical prediction
        props = self.predict_amp_properties(seq)
        result["properties"] = props

        # 2. Database similarity (if loaded)
        db_score = 0.0
        best_match: Optional[Dict[str, Any]] = None
        if self._amp_sequences and n >= 10:
            best_identity = 0.0
            best_id = ""
            # Simple identity — real version would use BLAST or k-mer Jaccard
            for amp_id, amp_seq in self._amp_sequences:
                if len(amp_seq) < n * 0.5 or len(amp_seq) > n * 2:
                    continue
                # K-mer Jaccard similarity (k=3)
                amp_kmers = set(
                    amp_seq[i:i+3] for i in range(len(amp_seq) - 2)
                )
                seq_kmers = set(seq[i:i+3] for i in range(n - 2))
                if not amp_kmers or not seq_kmers:
                    continue
                jaccard = len(amp_kmers & seq_kmers) / len(amp_kmers | seq_kmers)
                if jaccard > best_identity:
                    best_identity = jaccard
                    best_id = amp_id
            if best_identity > 0.2:
                db_score = best_identity
                best_match = {"id": best_id, "jaccard_similarity": round(best_identity, 4)}

        # 3. Combined score
        combined = max(props["score"], db_score)

        if best_match:
            result["best_db_match"] = best_match

        # AMP threshold
        is_amp = combined >= 0.5 and 5 <= n <= 150
        result["is_amp"] = is_amp
        result["confidence"] = round(combined, 4)

        return result


# ---------------------------------------------------------------------------
# Main classifier
# ---------------------------------------------------------------------------

class SmallPeptideClassifier:
    """Classify candidate sORF-encoded peptides into functional categories.

    Integrates:
    - SSP family classifier (motif matching)
    - Location classifier (GFF-based genomic context)
    - AMP classifier (physico-chemical + database)
    - miPEP detection (pri-miRNA scanning)
    - lncORF context detection

    Parameters
    ----------
    gff_path : str, optional
        Path to GFF/GTF annotation file for location-based classification.
    amp_db_path : str, optional
        Path to AMP database (FASTA format).
    custom_motifs : dict, optional
        Additional SSP family motif patterns.
    """

    def __init__(
        self,
        gff_path: Optional[str] = None,
        amp_db_path: Optional[str] = None,
        custom_motifs: Optional[Dict[str, str]] = None,
    ) -> None:
        self.ssp_classifier = SSPFamilyClassifier(custom_motifs)
        self.location_classifier = LocationClassifier(gff_path)
        self.amp_classifier = AMPClassifier(amp_db_path)

    # ---- Main classification entry point ----

    def classify_sorf(
        self,
        seq_id: str,
        sequence: str,
        gff_record: Optional[Dict[str, Any]] = None,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run full classification on a single sORF.

        Parameters
        ----------
        seq_id : str
            sORF identifier.
        sequence : str
            Amino acid sequence.
        gff_record : dict, optional
            Pre-parsed GFF record with 'feature_type', 'chrom', etc.
        chrom, start, end, strand : optional
            Genomic coordinates (alternative to gff_record).

        Returns
        -------
        dict with classification results.
        """
        seq = sequence.strip().upper()
        length = len(seq)

        result = ClassificationResult(
            seq_id=seq_id,
            sequence=seq,
            length_aa=length,
        )

        # ---- 1. SSP family motif matching ----
        ssp_result = self.ssp_classifier.classify(seq)
        if ssp_result["best_family"]:
            result.ssp_family = ssp_result["best_family"]
            result.categories[CATEGORY_SSP] = ssp_result["confidence"]
            result.notes.append(
                f"SSP motif match: {ssp_result['best_family']} "
                f"(confidence={ssp_result['confidence']:.3f})"
            )

        # ---- 2. Location context ----
        if gff_record:
            ctx = LocationClassifier.classify_by_context(gff_record)
        elif chrom is not None and start is not None and end is not None:
            ctx_result = self.location_classifier.classify_location(
                chrom=chrom, start=start, end=end, strand=strand or "+"
            )
            ctx = ctx_result["context"]
        else:
            ctx = "intergenic"

        result.location_context = ctx

        # Assign category based on location
        if ctx == "utr5":
            result.categories[CATEGORY_UORF] = 0.8
            result.notes.append("Located in 5' UTR → classified as uORF")
        elif ctx == "mirna":
            result.categories[CATEGORY_MIPEP] = 0.8
            result.notes.append("Located in miRNA precursor → classified as miPEP")
        elif ctx in ("lncrna",):
            result.categories[CATEGORY_LNCORF] = 0.8
            result.notes.append("Located in lncRNA → classified as lncORF")

        # ---- 3. AMP classification ----
        if MIN_LENGTH[CATEGORY_AMP] <= length <= MAX_LENGTH[CATEGORY_AMP]:
            amp_result = self.amp_classifier.classify(seq)
            if amp_result["is_amp"]:
                result.categories[CATEGORY_AMP] = amp_result["confidence"]
                result.amp_match = amp_result.get("best_db_match")
                result.notes.append(
                    f"AMP detected (confidence={amp_result['confidence']:.3f})"
                )

        # ---- 4. Default: if no category assigned, flag as intergenic sORF ----
        if not result.categories:
            if length >= 10:
                result.categories["novel_sORF"] = 0.5
                result.notes.append("Novel intergenic sORF (no known category)")

        return result.to_dict()

    # ---- Batch processing ----

    def classify_batch(
        self,
        records: List[Dict[str, Any]],
        n_jobs: int = 1,
    ) -> List[Dict[str, Any]]:
        """Classify multiple sORFs.

        Parameters
        ----------
        records : list of dict
            Each dict must have 'seq_id' and 'sequence'.
            May also have 'gff_record', 'chrom', 'start', 'end', 'strand'.
        n_jobs : int
            Number of workers (1 = serial).

        Returns
        -------
        list of dict
        """
        results: List[Dict[str, Any]] = []
        for rec in records:
            res = self.classify_sorf(
                seq_id=rec["seq_id"],
                sequence=rec["sequence"],
                gff_record=rec.get("gff_record"),
                chrom=rec.get("chrom"),
                start=rec.get("start"),
                end=rec.get("end"),
                strand=rec.get("strand"),
            )
            results.append(res)
        return results

    # ---- Report generation ----

    def generate_report(
        self,
        classifications: List[Dict[str, Any]],
        output_path: Optional[str] = None,
        format: str = "json",
    ) -> str:
        """Generate a summary report from classification results.

        Produces:
        - Per-category counts and proportions
        - SSP family breakdown
        - Location context distribution
        - Length distribution statistics

        Parameters
        ----------
        classifications : list of dict
            Results from classify_batch() or multiple classify_sorf() calls.
        output_path : str, optional
            Write report to this file.
        format : str
            'json' or 'tsv'.

        Returns
        -------
        str
            Report content (JSON string or TSV).
        """
        n_total = len(classifications)

        # Category counts
        cat_counts: Counter[str] = Counter()
        ssp_families: Counter[str] = Counter()
        locations: Counter[str] = Counter()
        lengths: List[int] = []
        confidences: Dict[str, List[float]] = defaultdict(list)

        for res in classifications:
            cats = res.get("categories", {})
            for cat, conf in cats.items():
                cat_counts[cat] += 1
                if cat == CATEGORY_SSP:
                    fam = res.get("ssp_family")
                    if fam:
                        ssp_families[fam] += 1
                confidences[cat].append(conf)

            loc = res.get("location_context", "unknown")
            locations[loc] += 1

            lengths.append(res.get("length_aa", 0))

        # Build report
        report: Dict[str, Any] = {
            "summary": {
                "total_sorfs_classified": n_total,
                "n_categories": len(cat_counts),
            },
            "category_counts": dict(cat_counts),
            "category_frequencies": {
                cat: round(count / n_total, 4) for cat, count in cat_counts.items()
            },
            "ssp_family_breakdown": dict(ssp_families),
            "location_distribution": dict(locations),
            "length_statistics": {
                "min": min(lengths) if lengths else 0,
                "max": max(lengths) if lengths else 0,
                "mean": round(np.mean(lengths), 2) if lengths else 0,
                "median": float(np.median(lengths)) if lengths else 0,
            },
            "category_confidence": {
                cat: {
                    "mean": round(np.mean(confs), 4),
                    "min": round(min(confs), 4),
                    "max": round(max(confs), 4),
                }
                for cat, confs in confidences.items()
            },
        }

        # Add multi-category stats
        n_multi = sum(
            1 for res in classifications
            if len(res.get("categories", {})) > 1
        )
        report["summary"]["n_multi_category"] = n_multi
        report["summary"]["fraction_multi_category"] = round(n_multi / n_total, 4) if n_total else 0

        if format == "tsv":
            return self._report_to_tsv(report, classifications)
        else:
            content = json.dumps(report, indent=2)
            if output_path:
                with open(output_path, "w") as f:
                    f.write(content)
                    f.write("\n")
                logger.info("Report written to %s", output_path)
            return content

    @staticmethod
    def _report_to_tsv(
        report: Dict[str, Any],
        classifications: List[Dict[str, Any]],
    ) -> str:
        """Generate TSV report with one row per sORF."""
        lines: List[str] = []
        # Header
        lines.append(
            "\t".join([
                "seq_id", "length_aa", "categories",
                "ssp_family", "location_context", "amp_match", "notes"
            ])
        )
        for res in classifications:
            cats = ";".join(
                f"{k}={v:.3f}" for k, v in res.get("categories", {}).items()
            )
            amp = json.dumps(res.get("amp_match") or "")
            notes = "; ".join(res.get("notes", []))
            lines.append(
                "\t".join([
                    res.get("seq_id", ""),
                    str(res.get("length_aa", "")),
                    cats,
                    res.get("ssp_family") or "",
                    res.get("location_context", ""),
                    amp,
                    notes,
                ])
            )
        return "\n".join(lines)

    # ---- SLURM support ----

    @staticmethod
    def chunk_for_slurm(
        records: List[Dict[str, Any]],
        n_chunks: int,
        output_dir: str,
        prefix: str = "classify_chunk",
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
    def merge_classifications(
        chunk_dir: str,
        output_path: str,
        pattern: str = "classify_chunk_*.json.results",
    ) -> str:
        """Merge per-chunk classification results into a single file."""
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
            "Merged classifications from %s → %s (%d records)",
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
        description="PeptSesame Layer3 — Small peptide classification"
    )
    parser.add_argument("--sequence", "-s", required=True, help="Amino acid sequence")
    parser.add_argument("--seq-id", "-i", default="test", help="Sequence identifier")
    parser.add_argument(
        "--context", "-c", default="intergenic",
        choices=["intergenic", "utr5", "mirna", "lncrna", "cds", "intron"],
        help="Genomic context",
    )
    args = parser.parse_args()

    classifier = SmallPeptideClassifier()

    result = classifier.classify_sorf(
        seq_id=args.seq_id,
        sequence=args.sequence,
        gff_record={"feature_type": args.context},
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
