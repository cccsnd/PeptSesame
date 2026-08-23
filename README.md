# PeptSesame

An evidence-layered, empirically calibrated prioritization framework for
sORF-encoded small-peptide candidates in plant genomes — from genome
sequence alone.

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

Small open reading frames (sORFs, 30–300 nt) encode functional micropeptides
that regulate plant growth, development, and stress responses. Their discovery
traditionally requires Ribo-seq or mass spectrometry data, which are
unavailable for most plant species. **PeptSesame** replaces experimental
evidence with a calibrated, rule-based computational scoring system, enabling
candidate prioritization for sORF-encoded peptides from any genome assembly.
The output is a **candidate-prioritization resource, not a catalogue of
confirmed peptides**: motif matches are reported as motif-compatible candidates
with transparent priority tiers, benchmarked recall/FPR, and external
mass-spectrometry-based coverage checks.

## Four-Layer Architecture

```
Layer 1: Six-frame translation + CDS-overlap filtering → sORF catalog
Layer 2: Multi-evidence scoring (5 rule channels) → priority score tiers
Layer 3: Motif-compatible classification (8 SSP families + IDA-like flag)
Layer 4: Functional annotation (GO / KEGG / cross-species conservation)
```

The five rule-based evidence channels are: sequence features (k-mer
composition, GC bias), cross-species conservation, RNA-seq expression,
structural features (signal peptide, transmembrane, cysteine richness), and
known motif matching. An optional machine-learning coding-potential channel is
defined but contributes a neutral value in the default (rule-based) mode.

## Repository layout

- `pipeline/` — importable Python package (installable via `pip install -e .`); provides the four-layer modules.
- `scripts/workflow/` — executable end-to-end run scripts that drive the package per species.

## Publication Pipeline

The reproducible production pipeline runs end-to-end per species:

```
scripts/workflow/run_layer1.py <species>     # six-frame translation → sORF BED/FASTA
scripts/workflow/run_layer2.py <species>     # multi-evidence scoring → scored_sorfs.tsv
scripts/workflow/run_layer3.py <species>     # motif classification → classified_sorfs.tsv
scripts/workflow/layer4_prep.py / layer4_run.py / layer4_extract.py  # SSP extraction + BLASTP prep/run
scripts/workflow/layer4_blastp.py <species>  # cross-species conservation
scripts/workflow/layer4_novel.py <species>   # novelty vs annotated proteins
scripts/workflow/layer4_core.py              # conserved × novel core set
scripts/workflow/run_expression.py           # sORF→gene mapping + DE (Yu11 charcoal rot)
scripts/workflow/benchmark_*.py              # motif recall/PPV, weight sensitivity, null model
scripts/workflow/tables*.py / generate_figures.py   # manuscript tables & figures
```

Input genomes/GFFs are listed in `config/pipeline_config.example.yaml`.
Internal paths are resolved via the `PEPTSESAME_ROOT` environment variable
(defaults to the repository root).

## Installation

```bash
pip install -e .          # core dependencies (biopython, numpy, pandas, scipy)
pip install -e ".[ml]"    # optional: LightGBM coding-potential channel
```

Python ≥ 3.10. BLAST+ (`blastp`) is required for the cross-species layer.

## Usage

```bash
export PEPTSESAME_ROOT=$(pwd)
mkdir -p data/genomes results_v2
# 1) place genome FASTA + GFF under data/ and edit config/pipeline_config.example.yaml
python scripts/workflow/run_layer1.py Yu11
python scripts/workflow/run_layer2.py Yu11
python scripts/workflow/run_layer3.py Yu11
```

## Citation

[To be added upon publication]

## License

MIT

## Key Results (18 plant genomes)

- Raw-ORF density under a fixed six-frame scan is a stable,
  composition-determined baseline (12,656–15,778 ORFs per Mb, CV = 6.0%),
  confirmed by nucleotide-composition-preserving null models
  (observed/random ratio 0.93–1.12)
- Motif-compatible SSP candidate resource: 3,020–6,922 candidates per sesame
  genome across eight conserved families (CLE, RALF, CEP, PSK, PSY1, IDA,
  EPFL, RGF), with 92.3% recall (12/13) against Arabidopsis gold-standard
  members (four families with curated gold standards), ≤0.06–0.11% FPR on
  random and protein-background peptides, and 73–90% coordinate-level
  coverage of independent MS-validated Arabidopsis peptides
- 30 SSP-motif-associated genes differentially expressed at BH FDR < 0.05 in
  stems during charcoal rot (*Macrophomina phaseolina*) infection (plus 221
  effect-size-selected nominal loci), enriched for conserved-novel
  candidates (OR = 3.08, P = 1.8×10⁻¹⁰; matched-background permutation
  P = 0.0044)

## Installation

```bash
# Conda environment
conda create -n pept python=3.11
conda activate pept
pip install -e .

# Dependencies
pip install biopython numpy pandas scikit-learn lightgbm pyyaml
```

## Quick Start

```bash
# Layer 1: six-frame translation
python scripts/run_layer1.py \
  --genome genome.fasta \
  --output results/sorfs.bed

# Layer 2: scoring
python scripts/run_layer2.py \
  --bed results/sorfs.bed \
  --fasta results/sorfs.fa \
  --output results/scored_sorfs.tsv

# Layer 3: classification
python scripts/run_layer3.py \
  --scored results/scored_sorfs.tsv \
  --peptide-fasta results/sorfs.fa \
  --output results/classified_sorfs.tsv \
  --genome-size 305

# Layer 4: functional annotation
python scripts/run_layer4.py \
  --classified results/classified_sorfs.tsv \
  --output results/functional \
  --genome-size 305
```

## SLURM Support

Batch processing for large genomes (100 Mb – 3 Gb) via SLURM array jobs:
```bash
sbatch scripts/slurm_batch_layer234.sh   # 11 species, array 0-10
```

## Performance

Rule-based scoring: ~11,000 sORFs/s per core. A 3-Gb genome (45 M sORFs)
completes all four layers in ~4 h on 2 cores.

## Citation

PeptSesame prioritizes motif-compatible small-peptide candidates in sesame
genomes through calibrated sequence and expression evidence (manuscript in
preparation).

## Code availability

The PeptSesame source code is publicly available at
[https://github.com/<org>/peptsesame](https://github.com/<org>/peptsesame)
under the MIT license, and archived at Zenodo
(https://doi.org/10.5281/zenodo.22064668).

- **Pipeline + scripts:** this repository (no genome data committed)
- **Full results (18-genome motif-compatible candidate tables, 30 FDR-DEGs /
  221 effect-size-selected loci, PTM propensities):**
  deposited at Figshare/Zenodo (see manuscript Data Availability)
- **Genome assemblies:** public databases (NCBI / NGDC / TAIR / IRGSP /
  MaizeGDB); accession lists in `config/pipeline_config.example.yaml`

For reproducibility, every pipeline stage in `scripts/workflow/` runs a
built-in self-check (row counts, score ranges, family-count reconciliation),
and every figure in the manuscript is generated by a script that reads
directly from the deposited data tables.

## License

MIT License — see [LICENSE](LICENSE).

## Contact

For questions or collaboration: chczhima@163.com; hengchuncao@gmail.com
