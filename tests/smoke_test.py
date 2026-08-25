#!/usr/bin/env python3
"""Smoke test for PeptSesame: exercises the Layer 1-4 entry points on a tiny demo genome.

Run:  python tests/smoke_test.py
Requires: biopython, numpy, pandas (see environment.yml)
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEMO_FASTA = """>chr1_demo
ATGGCTAAATGTACGGATCCTTAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGATCC
ATGGCTAAATGTACGGATCCTTAGGATCCATGGCTAAATGTACGGATCCTTAGGATCCATGGCTAAAT
GTACGGATCCTTAGGATCCATGGCTAAATGTACGGATCCTTAGGATCC
>chr2_demo
ATGGCTAAATGTACGGATCCTTAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGATCC
ATGGCTAAATGTACGGATCCTTAGGATCCATGGCTAAATGTACGGATCCTTAGGATCC
"""


def run(cmd, cwd):
    # ROOT 由脚本从自身位置推断 (仓库根); --outdir 控制输出位置, 不污染仓库
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"FAILED: {' '.join(cmd)}\n{r.stdout}\n{r.stderr}")
    return r


def main():
    checks = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "demo.fa").write_text(DEMO_FASTA)
        (td / "demo.gff").write_text(
            "chr1_demo\tdemo\tCDS\t1\t200\t.\t+\t.\tID=CDS_demo_001\n")

        # Layer 1
        run([sys.executable, str(ROOT / "scripts/workflow/run_layer1.py"),
             "demo", "demo.fa", "demo.gff", "--outdir", str(td)], td)
        bed = td / "results/01_layer1_sixframe/demo/sorfs.bed"
        if not bed.exists():
            raise RuntimeError("Layer 1 did not produce sorfs.bed")
        n = sum(1 for _ in open(bed)) - 1
        assert n > 0, "Layer 1 produced zero sORFs"
        checks += 1
        print(f"[ok] Layer 1: {n} sORFs")

        # Layer 2
        run([sys.executable, str(ROOT / "scripts/workflow/run_layer2.py"), "demo"], td)
        scored = td / "results/02_layer2_scoring/demo/scored_sorfs.tsv"
        assert scored.exists(), "Layer 2 did not produce scored_sorfs.tsv"
        checks += 1
        print("[ok] Layer 2: scored_sorfs.tsv")

        # Layer 3
        run([sys.executable, str(ROOT / "scripts/workflow/run_layer3.py"), "demo"], td)
        cls = td / "results/03_layer3_classify/demo/classified_sorfs.tsv"
        assert cls.exists(), "Layer 3 did not produce classified_sorfs.tsv"
        checks += 1
        print("[ok] Layer 3: classified_sorfs.tsv")

    # Layer 4 (cross-species conservation / novelty / shortlists) requires the
    # full 18-genome SSP library and is exercised by the real-data reruns, not
    # by this smoke test (see README "Reproducibility").
    print(f"\nSMOKE TEST PASSED ({checks} checks: Layer 1-3 core pipeline)")


if __name__ == "__main__":
    main()
