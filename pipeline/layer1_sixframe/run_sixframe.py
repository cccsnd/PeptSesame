#!/usr/bin/env python3
"""
run_sixframe.py — Standalone runner for PeptSesame Layer1 six-frame translation.

Can be used in three modes:

**Mode 1 — Real pipeline run** (specify your genome + GFF3)::

    python run_sixframe.py \\
        --genome /path/to/genome.fasta \\
        --gff /path/to/annotations.gff3 \\
        --out-dir <tmpdir>/pept_layer1_out \\
        --min-orf 30 --max-orf 300 \\
        --jobs 16

**Mode 2 — Demo with a synthetic mini-genome** (no input files needed)::

    python run_sixframe.py --demo

**Mode 3 — Verbose / debug**::

    python run_sixframe.py --genome ... --gff ... --verbose
"""

from __future__ import annotations

import argparse
import logging
import random
import sys
from pathlib import Path

# Ensure the parent of pipeline/ is on sys.path so the imports below work.
# When installed as a package this won't be needed, but for development
# we add the project root to the path.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]  # pipeline/../.. → pept/
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from pipeline.layer1_sixframe.sixframe import SixFrameTranslator  # noqa: E402

logger = logging.getLogger("run_sixframe")


# ---------------------------------------------------------------------------
# Synthetic mini-genome for demo mode
# ---------------------------------------------------------------------------

def _make_demo_genome(path: Path) -> None:
    """Write a tiny synthetic genome FASTA for demo/testing.

    Contains two contigs with a mix of intergenic sORFs and CDS regions
    so the overlap filter can be exercised.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # Contig 1 — 10 kbp with several embedded sORFs + one CDS-like region
    # We'll create a sequence that has known ORFs at predictable locations
    # Use deterministic random for reproducibility
    rng = random.Random(42)

    bases = ["A", "C", "G", "T"]
    contig1 = []

    # Region 0-1000: random intergenic
    contig1.append("".join(rng.choices(bases, k=1000)))

    # Region 1000-1120: an in-frame sORF (120 nt = 40 aa) with ATG start
    sorf1 = "ATG" + "".join(rng.choices(bases, k=114)) + "TGA"
    contig1.append(sorf1)

    # 1120-2000: random
    contig1.append("".join(rng.choices(bases, k=880)))

    # 2000-2600: a "known CDS" region (600 nt = 200 aa)
    cds1 = "ATG" + "".join(rng.choices(bases, k=594)) + "TGA"
    contig1.append(cds1)

    # 2600-3000: random
    contig1.append("".join(rng.choices(bases, k=400)))

    # 3000-3150: another sORF (150 nt = 50 aa)
    sorf2 = "ATG" + "".join(rng.choices(bases, k=144)) + "TAG"
    contig1.append(sorf2)

    # 3150-5000: random
    contig1.append("".join(rng.choices(bases, k=1850)))

    # 5000-5150: sORF overlapping the CDS region at 2000-2600? No, far away.
    # Let's put one right inside the CDS at 5200-5350
    contig1.append("".join(rng.choices(bases, k=200)))
    sorf3 = "ATG" + "".join(rng.choices(bases, k=144)) + "TAA"  # inside CDS range
    contig1.append(sorf3)

    # remainder
    contig1.append("".join(rng.choices(bases, k=4650)))  # total ~10k

    seq1 = "".join(contig1)[:10000]

    # Contig 2 — 5 kbp, simpler
    contig2 = "".join(rng.choices(bases, k=5000))

    with open(path, "w") as f:
        f.write(">chr1_demo\n")
        for i in range(0, len(seq1), 80):
            f.write(seq1[i : i + 80] + "\n")
        f.write(">chr2_demo\n")
        for i in range(0, len(contig2), 80):
            f.write(contig2[i : i + 80] + "\n")

    logger.info("Demo genome written to %s", path)


def _make_demo_gff(path: Path) -> None:
    """Write a demo GFF3 with a single CDS feature for exercise."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(
            "chr1_demo\tdemo\tCDS\t2001\t2600\t.\t+\t.\t"
            "ID=CDS_demo_001;Name=demo_gene\n"
        )
    logger.info("Demo GFF3 written to %s", path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    p = argparse.ArgumentParser(
        description="PeptSesame Layer 1 — Six-frame translation & sORF filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input
    p.add_argument("--genome", type=str, help="Path to genome FASTA file")
    p.add_argument("--gff", type=str, default=None, help="Path to CDS GFF3 file")

    # Output
    p.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Output directory (default: current working directory)",
    )
    p.add_argument("--out-bed", type=str, default=None, help="Output BED filename")
    p.add_argument(
        "--out-fasta", type=str, default=None, help="Output peptide FASTA filename"
    )

    # Parameters
    p.add_argument(
        "--min-orf",
        type=int,
        default=30,
        help="Minimum ORF length in nt (default: 30 = 10 aa)",
    )
    p.add_argument(
        "--max-orf",
        type=int,
        default=300,
        help="Maximum ORF length in nt (default: 300 = 100 aa)",
    )
    p.add_argument(
        "--strand-agnostic",
        action="store_true",
        help="Filter overlapping CDS regardless of strand "
        "(default: only same-strand CDS overlap is filtered)",
    )
    p.add_argument(
        "--keep-partial",
        action="store_true",
        help="Include ORFs that run to contig end without a stop codon",
    )

    # Performance
    p.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Number of parallel workers (default: all CPU cores)",
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Chunk size (bp) for processing large contigs",
    )

    # Temp
    p.add_argument(
        "--tmp-dir",
        type=str,
        default=tempfile.gettempdir(),
        help="Temporary directory (default: system temp directory)",
    )

    # Mode
    p.add_argument(
        "--demo",
        action="store_true",
        help="Run with a synthetic mini-genome for testing",
    )

    # Logging
    p.add_argument(
        "-v", "--verbose", action="store_true", help="Increase log verbosity"
    )

    return p


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    # Logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    # --- Demo mode ---
    if args.demo:
        tmp_dir = Path(args.tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        demo_fasta = tmp_dir / "demo_genome.fasta"
        demo_gff = tmp_dir / "demo_annotations.gff3"
        _make_demo_genome(demo_fasta)
        _make_demo_gff(demo_gff)
        args.genome = str(demo_fasta)
        args.gff = str(demo_gff)
        if args.out_dir is None:
            args.out_dir = str(tmp_dir / "demo_output")
        logger.info("=== DEMO MODE ===")
        logger.info("Genome: %s", demo_fasta)
        logger.info("GFF:    %s", demo_gff)
        logger.info("Output: %s", args.out_dir)

    # --- Validate ---
    if not args.genome:
        parser.error("Either --genome or --demo is required")
    genome_path = Path(args.genome)
    if not genome_path.exists():
        parser.error(f"Genome FASTA not found: {genome_path}")

    gff_path = Path(args.gff) if args.gff else None
    if gff_path is not None and not gff_path.exists():
        parser.error(f"GFF3 not found: {gff_path}")

    # --- Output paths ---
    out_dir = Path(args.out_dir) if args.out_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_bed = out_dir / (args.out_bed or "sorfs.bed")
    out_fasta = out_dir / (args.out_fasta or "sorfs_peptide.fasta")

    # --- Build translator ---
    st = SixFrameTranslator(
        genome_fasta=genome_path,
        cds_gff=gff_path,
        out_bed=out_bed,
        out_pep_fasta=out_fasta,
        min_orf_len=args.min_orf,
        max_orf_len=args.max_orf,
        strand_aware_overlap=not args.strand_agnostic,
        n_jobs=args.jobs,
        tmp_dir=args.tmp_dir,
        chunk_size=args.chunk_size,
        keep_partial_orfs=args.keep_partial,
    )

    # --- Run ---
    try:
        orfs = st.run()
    except Exception:
        logger.exception("Pipeline failed")
        sys.exit(1)

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"  PeptSesame Layer 1 — Complete")
    print(f"  Total sORF candidates: {len(orfs)}")
    print(f"  Output BED  : {out_bed}")
    print(f"  Output FASTA: {out_fasta}")
    print(f"{'='*60}\n")

    # Quick stats
    if orfs:
        strands = {}
        frames = {}
        lengths = [o.length_nt for o in orfs]
        for o in orfs:
            strands[o.strand] = strands.get(o.strand, 0) + 1
            frames[o.frame] = frames.get(o.frame, 0) + 1
        print("  Strand distribution:")
        for s in ("+", "-"):
            print(f"    {s}: {strands.get(s, 0)}")
        print("  Frame distribution:")
        for f in (0, 1, 2):
            print(f"    Frame {f}: {frames.get(f, 0)}")
        print(f"  Length range: {min(lengths)}–{max(lengths)} nt")
        print(f"  Mean length: {sum(lengths)/len(lengths):.1f} nt")


if __name__ == "__main__":
    main()
