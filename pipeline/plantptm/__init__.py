"""PlantPTM web batch client + result parsing (Paper-2 Evidence-D module).

PlantPTM (Molecular Plant 2026, 9 PTM types, 6 training species) — batch
post-translational modification site prediction for plant peptides.

- client: CSRF-authenticated batch submission + task polling
- parse:  raw batch JSON -> site rows, protein summary, sequon sanity

Protocol & pitfalls (clean FASTA IDs only; urlencoded POST; full-length
submission; confidence-stratified sequon validation):
see skill plant-small-peptide-mining -> references/plantptm-batch-submission.md
"""
from __future__ import annotations

from .client import (
    BASE,
    PTM_TYPES,
    THRESHOLDS,
    UA,
    PlantPTMClient,
    read_fasta,
)
from .parse import (
    PTM_CN,
    PTM_RES,
    load_raw_rows,
    sanity_ngly_sequon,
)

__all__ = [
    "BASE", "PTM_TYPES", "THRESHOLDS", "UA", "PlantPTMClient", "read_fasta",
    "PTM_CN", "PTM_RES", "load_raw_rows", "sanity_ngly_sequon",
]
