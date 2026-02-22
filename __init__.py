"""
LIR Motif Finder
================
Find LC3-Interacting Region (LIR) motifs in intrinsically disordered
regions predicted by metapredict.

Quick start
-----------
>>> from lirmotiffinder import run_analysis
>>> results = run_analysis("proteins.fasta")
>>> for result in results:
...     for hit in result.motif_hits:
...         print(hit.protein_id, hit.motif_sequence, hit.motif_start_protein)
"""

from .core import (
    run_analysis,
    analyze_protein,
    predict_disorder,
    find_lir_in_region,
    DisorderedRegion,
    MotifHit,
    ProteinResult,
)

__all__ = [
    "run_analysis",
    "analyze_protein",
    "predict_disorder",
    "find_lir_in_region",
    "DisorderedRegion",
    "MotifHit",
    "ProteinResult",
]

__version__ = "0.1.0"
__author__ = "Your Name"
