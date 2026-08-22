"""
sixframe.py — Six-frame genome translation and sORF filtering.

Provides :class:`SixFrameTranslator` which:
1. Loads a genome FASTA and known CDS annotations (GFF3).
2. Performs six-frame translation (forward 3 frames + reverse 3 frames).
3. Extracts ORFs between ``min_length`` and ``max_length`` nucleotides.
4. Filters out ORFs overlapping known CDS regions.
5. Exports a BED file and a peptide FASTA file.

Typical command-line usage is provided by the sibling script ``run_sixframe.py``.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import tempfile
import warnings
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings

# Suppress the Biopython partial-codon warning — this is expected in
# six-frame translation where the first frame offset often leaves trailing
# partial codons.  We filter by message text so only the harmless one is
# silenced.
warnings.filterwarnings(
    "ignore",
    message="Partial codon",
    category=Warning,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default minimum ORF length in nucleotides (10 aa × 3 nt).
DEFAULT_MIN_ORF_LEN = 30
#: Default maximum ORF length in nucleotides (100 aa × 3 nt).
DEFAULT_MAX_ORF_LEN = 300
#: Codon table to use for translation (1 = standard / universal).
CODON_TABLE_ID = 1
#: Pattern for extracting a coding ORF from a translated AA string: start-M
#: followed by any non-stop residues, ending at stop (*) or contig boundary.
#: Uses lookbehind so overlapping ORFs (rare in short peptides) are still
#: observable; in practice we scan frame-by-frame so greedy is fine.
_ORF_RE = re.compile(r"M[^*]*\*")
#: Number of frames per strand.
_FRAMES = (0, 1, 2)
#: BED header line.
_BED_HEADER = "#chrom\tstart\tend\tname\tscore\tstrand\tframe\tlength_nt\tlength_aa\tpeptide_seq"

# ---------------------------------------------------------------------------
# Simple interval container
# ---------------------------------------------------------------------------


class _Interval:
    """A half-open interval [start, end) on a strand."""

    __slots__ = ("start", "end", "strand")

    def __init__(self, start: int, end: int, strand: str) -> None:
        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self) -> str:
        return f"_Interval({self.start}, {self.end}, {self.strand!r})"


# ---------------------------------------------------------------------------
# ORF record
# ---------------------------------------------------------------------------


class ORFRecord:
    """Data container for a single small ORF candidate.

    Parameters
    ----------
    contig : str
        Chromosome / scaffold name.
    start : int
        0-based genomic start coordinate (inclusive).
    end : int
        0-based genomic end coordinate (exclusive).
    strand : str
        ``'+'`` or ``'-'``.
    frame : int
        Frame on that strand (0, 1, or 2).
    peptide_seq : str
        Translated amino-acid sequence.
    length_nt : int
        ORF length in nucleotides (``end - start``).
    """

    __slots__ = (
        "contig",
        "start",
        "end",
        "strand",
        "frame",
        "peptide_seq",
        "length_nt",
    )

    def __init__(
        self,
        contig: str,
        start: int,
        end: int,
        strand: str,
        frame: int,
        peptide_seq: str,
    ) -> None:
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.frame = frame
        self.peptide_seq = peptide_seq
        self.length_nt = end - start

    @property
    def length_aa(self) -> int:
        """Length of the translated peptide."""
        return len(self.peptide_seq)

    def to_bed_fields(self) -> Tuple[str, int, int, str, int, str, int, int, int, str]:
        """Return fields for the extended BED output.

        Columns: chrom, start, end, name, score, strand, frame, length_nt,
        length_aa, peptide_seq.
        """
        name = f"sORF_{self.contig}_{self.start}_{self.end}_{self.strand}_f{self.frame}"
        score = 0  # placeholder; can be replaced by Layer 2 scoring
        return (
            self.contig,
            self.start,
            self.end,
            name,
            score,
            self.strand,
            self.frame,
            self.length_nt,
            self.length_aa,
            self.peptide_seq,
        )

    def to_seq_record(self) -> SeqRecord:
        """Return a BioPython SeqRecord for FASTA export."""
        desc = (
            f"strand={self.strand} frame={self.frame} "
            f"length_nt={self.length_nt} length_aa={self.length_aa}"
        )
        return SeqRecord(
            Seq(self.peptide_seq),
            id=f"{self.contig}_{self.start}_{self.end}_{self.strand}_f{self.frame}",
            description=desc,
        )


# ---------------------------------------------------------------------------
# Interval-tree helper for overlap queries
# ---------------------------------------------------------------------------


class _IntervalTree:
    """Minimal interval tree for half-open genomic intervals.

    Supports add and overlap query.  Implemented as a sorted list with
    bisect-based lookups — sufficient for the moderate number of CDS
    features per contig in a typical plant genome (tens of thousands).
    """

    __slots__ = ("_starts", "_ends")

    def __init__(self) -> None:
        self._starts: List[int] = []
        self._ends: List[int] = []

    def add(self, start: int, end: int) -> None:
        """Add an interval.  Runs in O(1) amortised."""
        self._starts.append(start)
        self._ends.append(end)

    def finalize(self) -> None:
        """Sort by start coordinate (call once after all adds)."""
        if not self._starts:
            return
        pairs = sorted(zip(self._starts, self._ends))
        self._starts, self._ends = zip(*pairs) if pairs else ([], [])
        # Convert back to lists after zip produces tuples
        self._starts = list(self._starts)
        self._ends = list(self._ends)

    def overlaps(self, start: int, end: int) -> bool:
        """Return True if [start, end) overlaps any stored interval.

        Uses binary search (O(log n)).
        """
        import bisect

        if not self._starts:
            return False
        # Find the first interval whose start <= end-1 (i.e. could overlap)
        i = bisect.bisect_right(self._starts, end - 1) - 1
        if i < 0:
            return False
        # Check the intervals that could overlap backward
        while i >= 0 and self._ends[i] > start:
            if self._ends[i] > start:
                return True
            i -= 1
        return False


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------


class SixFrameTranslator:
    """Six-frame genome translation and sORF candidate filtering.

    Parameters
    ----------
    genome_fasta : str or Path
        Path to the genome FASTA file.
    cds_gff : str or Path or None
        Path to a GFF3 file containing CDS annotations for overlap filtering.
        Pass ``None`` to skip filter (not recommended).
    out_bed : str or Path or None
        Path for the output BED file.  Defaults to ``sorfs.bed`` in the
        current directory.
    out_pep_fasta : str or Path or None
        Path for the output peptide FASTA file.  Defaults to
        ``sorfs_peptide.fasta`` in the current directory.
    min_orf_len : int
        Minimum ORF length in nucleotides (default 30).
    max_orf_len : int
        Maximum ORF length in nucleotides (default 300).
    strand_aware_overlap : bool
        If ``True``, only filter sORFs that overlap a CDS on the **same**
        strand.  If ``False``, any overlapping CDS (any strand) triggers
        filtering.
    n_jobs : int
        Number of parallel worker processes.  Default: all CPU cores.
    tmp_dir : str or Path
        Directory for temporary files.  Defaults to the system temp directory.
        (fallback: system tempdir).
    chunk_size : int or None
        Bases per chunk when processing huge contigs (e.g., 1 Gbp).
        ``None`` = no chunking (process each contig as a whole).
    keep_partial_orfs : bool
        If True, include ORFs that extend to the end of a contig (no stop
        codon).  Such ORFs may be real but should be treated with caution.
    """

    def __init__(
        self,
        genome_fasta: str | Path,
        cds_gff: str | Path | None = None,
        out_bed: str | Path | None = None,
        out_pep_fasta: str | Path | None = None,
        min_orf_len: int = DEFAULT_MIN_ORF_LEN,
        max_orf_len: int = DEFAULT_MAX_ORF_LEN,
        strand_aware_overlap: bool = True,
        n_jobs: int | None = None,
        tmp_dir: str | Path = tempfile.gettempdir(),
        chunk_size: int | None = None,
        keep_partial_orfs: bool = False,
    ) -> None:
        self.genome_fasta = Path(genome_fasta)
        self.cds_gff = Path(cds_gff) if cds_gff else None
        self.out_bed = Path(out_bed or "sorfs.bed")
        self.out_pep_fasta = Path(out_pep_fasta or "sorfs_peptide.fasta")
        self.min_orf_len = min_orf_len
        self.max_orf_len = max_orf_len
        self.strand_aware_overlap = strand_aware_overlap
        self.n_jobs = n_jobs or cpu_count()
        self.tmp_dir = Path(tmp_dir)
        self.chunk_size = chunk_size
        self.keep_partial_orfs = keep_partial_orfs

        # Internal state — populated during run()
        self._genome: Dict[str, str] = {}
        self._cds_tree: Dict[str, _IntervalTree] = {}
        self._orfs: List[ORFRecord] = []

    # ---- Public API -------------------------------------------------------

    def run(self) -> List[ORFRecord]:
        """Execute the full pipeline.

        Returns the list of filtered :class:`ORFRecord` objects.
        """
        logger.info("PeptSesame Layer 1 — Six-frame translation")
        logger.info("  Genome FASTA : %s", self.genome_fasta)
        logger.info("  CDS GFF      : %s", self.cds_gff)
        logger.info("  Output BED   : %s", self.out_bed)
        logger.info("  Output FASTA : %s", self.out_pep_fasta)
        logger.info("  Jobs         : %d", self.n_jobs)
        logger.info("  ORF length   : %d–%d nt", self.min_orf_len, self.max_orf_len)

        # Step 1 — load genome
        self._load_genome()
        logger.info("Loaded %d contigs from genome FASTA", len(self._genome))

        # Step 2 — load CDS intervals (if provided)
        if self.cds_gff is not None and self.cds_gff.exists():
            self._load_cds_intervals()
            logger.info(
                "Loaded CDS intervals for %d contigs", len(self._cds_tree)
            )
        else:
            logger.info("No CDS annotation provided — skipping overlap filter")

        # Step 3 — six-frame translate and extract ORFs
        self._translate_and_extract()

        # Step 4 — filter by length (already done during extraction, but
        # we double-check here for safety)
        self._orfs = [o for o in self._orfs if self.min_orf_len <= o.length_nt <= self.max_orf_len]
        logger.info("After length filter: %d ORFs", len(self._orfs))

        # Step 5 — filter overlapping CDS
        if self._cds_tree:
            n_before = len(self._orfs)
            self._filter_overlapping_cds()
            n_removed = n_before - len(self._orfs)
            logger.info("Overlap filter removed %d ORFs (remaining: %d)", n_removed, len(self._orfs))
        else:
            logger.info("Skipping CDS overlap filter (no annotation)")

        # Step 6 — export
        self._write_bed()
        self._write_peptide_fasta()

        logger.info("Done. %d sORF candidates written.", len(self._orfs))
        return self._orfs

    # ---- Genome loading ---------------------------------------------------

    def _load_genome(self) -> None:
        """Load the genome FASTA into ``self._genome`` (contig → uppercase DNA)."""
        self._genome = {}
        for record in SeqIO.parse(str(self.genome_fasta), "fasta"):
            seq = str(record.seq).upper().replace("U", "T")
            # Clean — only keep standard IUPAC DNA bases
            self._genome[record.id] = seq
        if not self._genome:
            raise ValueError(f"No sequences found in {self.genome_fasta}")

    # ---- CDS interval loading ---------------------------------------------

    def _load_cds_intervals(self) -> None:
        """Parse GFF3 CDS features into per-contig interval trees.

        Handles the standard GFF3 column format:
            chr  source  CDS  start  end  score  strand  …
        Coordinates in GFF3 are 1-based inclusive; we convert to 0-based
        half-open for internal use.
        """
        cds_intervals: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        with open(self.cds_gff) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith(">"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature_type = parts[2].upper()
                if feature_type != "CDS":
                    continue
                chrom = parts[0]
                gff_start = int(parts[3])
                gff_end = int(parts[4])
                strand = parts[6]
                if strand not in ("+", "-", "."):
                    strand = "."
                # GFF3: 1-based inclusive → 0-based half-open
                start_0 = gff_start - 1
                end_0 = gff_end
                cds_intervals[chrom].append((start_0, end_0, strand))

        # Build interval trees
        self._cds_tree = {}
        for chrom, intervals in cds_intervals.items():
            tree = _IntervalTree()
            if self.strand_aware_overlap:
                # Only consider same-strand CDS for overlap
                # BUT: we store all and filter by strand at query time
                for s, e, _ in intervals:
                    tree.add(s, e)
            else:
                for s, e, _ in intervals:
                    tree.add(s, e)
            tree.finalize()
            self._cds_tree[chrom] = tree

        # Also build a strand-aware dict for the strand-aware overlap case
        if self.strand_aware_overlap:
            self._cds_tree_plus: Dict[str, _IntervalTree] = {}
            self._cds_tree_minus: Dict[str, _IntervalTree] = {}
            for chrom, intervals in cds_intervals.items():
                t_plus = _IntervalTree()
                t_minus = _IntervalTree()
                for s, e, strand in intervals:
                    if strand == "+":
                        t_plus.add(s, e)
                    elif strand == "-":
                        t_minus.add(s, e)
                t_plus.finalize()
                t_minus.finalize()
                self._cds_tree_plus[chrom] = t_plus
                self._cds_tree_minus[chrom] = t_minus

    # ---- Translation & ORF extraction (single-chunk worker) ----------------

    @staticmethod
    def _translate_contig(
        contig_name: str,
        dna_seq: str,
        min_len: int,
        max_len: int,
        keep_partial: bool,
    ) -> List[ORFRecord]:
        """Six-frame translate one contig and extract ORFs.

        Parameters
        ----------
        contig_name : str
            Contig / chromosome identifier.
        dna_seq : str
            Upper-case DNA sequence.
        min_len : int
            Minimum ORF length in **nucleotides**.
        max_len : int
            Maximum ORF length in **nucleotides**.
        keep_partial : bool
            Include ORFs lacking a stop codon (contig-end truncation).

        Returns
        -------
        list of ORFRecord
            ORF candidates found in this contig.
        """
        results: List[ORFRecord] = []

        # --- Forward strand -------------------------------------------------
        dna_seq_len = len(dna_seq)
        for frame in _FRAMES:
            # Translate the reading frame
            translated = str(
                Seq(dna_seq[frame:]).translate(table=CODON_TABLE_ID, to_stop=False)
            )
            # Find all ORFs starting with M and ending with *
            for match in _ORF_RE.finditer(translated):
                aa_start = match.start()
                aa_end = match.end()
                # Convert AA coordinates back to genomic nucleotide coordinates
                nt_start = frame + aa_start * 3
                # aa_end points *after* the stop codon in the AA string
                # The stop codon is included in the match, so:
                nt_end = frame + aa_end * 3  # inclusive of stop codon
                nt_len = nt_end - nt_start

                if nt_len < min_len or nt_len > max_len:
                    continue

                peptide = translated[aa_start : aa_end - 1]  # exclude stop
                if len(peptide) < 1:
                    continue

                results.append(
                    ORFRecord(
                        contig=contig_name,
                        start=nt_start,
                        end=nt_end,
                        strand="+",
                        frame=frame,
                        peptide_seq=peptide,
                    )
                )

            # Partial ORFs at contig end (no stop codon)
            if keep_partial:
                # Scan forward from each M until contig end
                i = 0
                while i < len(translated):
                    if translated[i] == "M":
                        # Find next stop or end
                        stop_idx = translated.find("*", i)
                        if stop_idx == -1:
                            # No stop → partial ORF running to contig end
                            aa_seq = translated[i:]
                            nt_start = frame + i * 3
                            nt_end = dna_seq_len
                            nt_len = nt_end - nt_start
                            if min_len <= nt_len <= max_len and len(aa_seq) >= 1:
                                results.append(
                                    ORFRecord(
                                        contig=contig_name,
                                        start=nt_start,
                                        end=nt_end,
                                        strand="+",
                                        frame=frame,
                                        peptide_seq=aa_seq,
                                    )
                                )
                            # No more M's after this point that could yield
                            # another partial, so break
                            break
                        else:
                            i = stop_idx + 1
                            continue
                    i += 1

        # --- Reverse strand -------------------------------------------------
        rev_seq = str(Seq(dna_seq).reverse_complement())
        rev_seq_len = len(rev_seq)
        for frame in _FRAMES:
            translated = str(
                Seq(rev_seq[frame:]).translate(table=CODON_TABLE_ID, to_stop=False)
            )
            for match in _ORF_RE.finditer(translated):
                aa_start = match.start()
                aa_end = match.end()
                nt_start_rev = frame + aa_start * 3
                nt_end_rev = frame + aa_end * 3
                nt_len = nt_end_rev - nt_start_rev

                if nt_len < min_len or nt_len > max_len:
                    continue

                peptide = translated[aa_start : aa_end - 1]
                if len(peptide) < 1:
                    continue

                # Convert reverse-strand coordinates back to forward-strand
                # genomic coordinates
                # The reverse-complement sequence runs 3'→5' on the forward
                # strand.  Position *i* in rev_seq corresponds to position
                # (dna_seq_len - 1 - i) on the forward strand.
                # We need the ORF on the FORWARD coordinate, end > start.
                end_fwd = dna_seq_len - nt_start_rev
                start_fwd = dna_seq_len - nt_end_rev

                results.append(
                    ORFRecord(
                        contig=contig_name,
                        start=start_fwd,
                        end=end_fwd,
                        strand="-",
                        frame=frame,
                        peptide_seq=peptide,
                    )
                )

            # Partial ORFs on reverse strand
            if keep_partial:
                i = 0
                while i < len(translated):
                    if translated[i] == "M":
                        stop_idx = translated.find("*", i)
                        if stop_idx == -1:
                            aa_seq = translated[i:]
                            nt_start_rev = frame + i * 3
                            nt_end_rev = rev_seq_len
                            nt_len = nt_end_rev - nt_start_rev
                            if min_len <= nt_len <= max_len and len(aa_seq) >= 1:
                                end_fwd = dna_seq_len - nt_start_rev
                                start_fwd = dna_seq_len - nt_end_rev
                                results.append(
                                    ORFRecord(
                                        contig=contig_name,
                                        start=start_fwd,
                                        end=end_fwd,
                                        strand="-",
                                        frame=frame,
                                        peptide_seq=aa_seq,
                                    )
                                )
                            break
                        else:
                            i = stop_idx + 1
                            continue
                    i += 1

        return results

    # ---- Parallel dispatching ---------------------------------------------

    def _translate_and_extract(self) -> None:
        """Dispatch contig translation across worker processes."""
        if self.n_jobs <= 1 or len(self._genome) == 1:
            # Sequential path (also used as the building block for chunking)
            self._orfs = []
            for contig_name, dna_seq in self._genome.items():
                logger.debug("Translating %s …", contig_name)
                self._orfs.extend(
                    self._translate_contig(
                        contig_name,
                        dna_seq,
                        self.min_orf_len,
                        self.max_orf_len,
                        self.keep_partial_orfs,
                    )
                )
            return

        # Parallel path
        args = [
            (
                contig_name,
                dna_seq,
                self.min_orf_len,
                self.max_orf_len,
                self.keep_partial_orfs,
            )
            for contig_name, dna_seq in self._genome.items()
        ]

        logger.info(
            "Translating %d contigs with %d workers …",
            len(args),
            self.n_jobs,
        )

        with Pool(self.n_jobs) as pool:
            results = pool.starmap(self._translate_contig, args)

        self._orfs = []
        for orf_list in results:
            self._orfs.extend(orf_list)

        logger.info("Translation complete: %d raw ORFs extracted", len(self._orfs))

    # ---- CDS overlap filter -----------------------------------------------

    def _filter_overlapping_cds(self) -> None:
        """Remove ORFs that overlap known CDS intervals."""
        filtered: List[ORFRecord] = []
        for orf in self._orfs:
            tree = self._cds_tree.get(orf.contig)
            if tree is None:
                filtered.append(orf)
                continue

            # Determine if we overlap
            if self.strand_aware_overlap:
                # Check specifically same-strand CDS
                if orf.strand == "+":
                    target_tree = self._cds_tree_plus.get(orf.contig)
                else:
                    target_tree = self._cds_tree_minus.get(orf.contig)

                if target_tree is None or not target_tree.overlaps(orf.start, orf.end):
                    filtered.append(orf)
            else:
                if not tree.overlaps(orf.start, orf.end):
                    filtered.append(orf)

        self._orfs = filtered

    # ---- Output writers ----------------------------------------------------

    def _write_bed(self) -> None:
        """Write the sORF candidates to an extended BED file."""
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        # Write to temp first, then rename for atomicity
        tmp = tempfile.NamedTemporaryFile(
            mode="w",
            dir=str(self.tmp_dir),
            suffix=".bed",
            delete=False,
        )
        tmp_path = tmp.name
        try:
            tmp.write(_BED_HEADER + "\n")
            for orf in sorted(self._orfs, key=lambda o: (o.contig, o.start)):
                fields = orf.to_bed_fields()
                tmp.write(
                    f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t"
                    f"{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t"
                    f"{fields[8]}\t{fields[9]}\n"
                )
            tmp.flush()
            os.fsync(tmp.fileno())
        finally:
            tmp.close()
        # Atomic rename to final destination
        shutil.move(tmp_path, str(self.out_bed))
        logger.info("Wrote %s (%d entries)", self.out_bed, len(self._orfs))

    def _write_peptide_fasta(self) -> None:
        """Write the peptide sequences to FASTA."""
        records = [orf.to_seq_record() for orf in self._orfs]
        self.out_pep_fasta.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records, str(self.out_pep_fasta), "fasta")
        logger.info("Wrote %s (%d sequences)", self.out_pep_fasta, len(records))
