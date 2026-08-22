# Contributing to PeptSesame

Thanks for your interest in PeptSesame! This document covers the conventions
we use so that issues and pull requests can be reviewed efficiently.

## Development setup

```bash
git clone https://github.com/<org>/peptsesame.git
cd peptsesame
conda create -n pept python=3.11
conda activate pept
pip install -e ".[ml]"

# Minimal smoke test
peptsesame --version
python scripts/self_check.py
```

## Code layout

```
pipeline/
  layer1_sixframe/   six-frame translation + CDS-overlap filtering
  layer2_scoring/    multi-evidence scoring (sequence/ML/conservation/expression/structure/motif)
  layer3_classify/   SSP-family / miPEP / AMP / uORF classification
  layer4_function/   GO / KEGG / stress-association annotation
  plantptm/          PlantPTM submission client + result parsing (post-translational modification layer)
scripts/             pipeline entry points, figure/table generation, analysis scripts
config/              YAML configuration (use pipeline_config.example.yaml as template)
```

## Conventions

1. **No hard-coded absolute paths in library code.** `pipeline/` modules must be
   path-agnostic; scripts read paths from `config/pipeline_config.yaml` or CLI
   arguments. Scripts that still need cluster paths must keep them behind a
   config lookup or clearly marked placeholder (see `run_cross_species.py`).
2. **Data files are not committed.** Genomes, results, and intermediate tables
   are ignored via `.gitignore`; only small example inputs are allowed under
   `examples/`.
3. **Figures are reproducible.** Every figure has a generating script that
   reads from TSV/CSV source files; do not hard-code counts inside plotting
   scripts (values go stale silently — see `generate_figures.py` for the
   current pattern).
4. **Self-check before pushing.** Run `python scripts/self_check.py` and make
   sure the pipeline imports cleanly (`python -c "import pipeline"`).

## Reporting issues

Please include:
- Python version and OS
- Exact command used
- Full traceback
- Input file format summary (2-3 lines from `head`)

## Pull requests

- Branch from `main`; one logical change per PR.
- Add a one-line summary in the PR description mapping to the affected layer.
- Run the smoke test above before requesting review.

## License

MIT — see [LICENSE](LICENSE). By contributing you agree to license your
contribution under the same terms.
