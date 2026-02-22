"""
Core logic for LIR motif detection in intrinsically disordered regions.

Motif definitions (based on canonical LIR/AIM motif literature):
  - Basic LIR:          [WFY]xx[LIV]
  - Extended LIR:       [DE][WFY]xx[LIV]
  - Acidic-extended:    [DE]{2}[WFY]xx[LIV]

Disordered regions are predicted using metapredict V2.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, asdict
from typing import Iterator

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False
    SeqIO = None  # type: ignore[assignment]
    SeqRecord = None  # type: ignore[assignment]

try:
    import metapredict as meta
except ImportError as exc:
    raise ImportError(
        "metapredict is required. Install it with:\n"
        "  pip install metapredict"
    ) from exc


# ---------------------------------------------------------------------------
# Built-in FASTA parser (fallback when BioPython is not available)
# ---------------------------------------------------------------------------
class _SimpleRecord:
    """Minimal stand-in for Bio.SeqRecord.SeqRecord."""
    __slots__ = ("id", "description", "_seq")

    def __init__(self, header: str, sequence: str) -> None:
        parts = header.split(None, 1)
        self.id = parts[0]
        self.description = header
        self._seq = sequence

    def __str__(self) -> str:
        return self._seq


def _parse_fasta(path: str):
    """Yield _SimpleRecord objects from a FASTA file."""
    header = ""
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield _SimpleRecord(header, "".join(chunks))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header:
        yield _SimpleRecord(header, "".join(chunks))


# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------
_MOTIF_PATTERNS: dict[str, re.Pattern[str]] = {
    "basic":            re.compile(r"(?=([WFY][A-Z]{2}[LIV]))"),
    "extended":         re.compile(r"(?=([DE][WFY][A-Z]{2}[LIV]))"),
    "acidic_extended":  re.compile(r"(?=([DE]{2}[WFY][A-Z]{2}[LIV]))"),
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class DisorderedRegion:
    """A contiguous disordered region within a protein."""
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive
    sequence: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class MotifHit:
    """A single LIR motif found within a disordered region."""
    protein_id: str
    protein_description: str
    protein_length: int
    motif_type: str
    motif_sequence: str
    # Positions relative to full protein (1-based for biologists)
    motif_start_protein: int
    motif_end_protein: int
    # Disordered region boundaries (1-based)
    disorder_start_protein: int
    disorder_end_protein: int
    disorder_length: int
    # Raw disorder score at motif position
    mean_disorder_score: float

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ProteinResult:
    """All results for a single protein."""
    protein_id: str
    protein_description: str
    protein_length: int
    n_disordered_regions: int
    n_disordered_residues: int
    motif_hits: list[MotifHit] = field(default_factory=list)

    @property
    def has_hits(self) -> bool:
        return len(self.motif_hits) > 0


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------
def predict_disorder(sequence: str, threshold: float = 0.5) -> list[DisorderedRegion]:
    """
    Run metapredict on *sequence* and return contiguous disordered regions.

    Parameters
    ----------
    sequence : str
        Amino-acid sequence (single-letter codes, uppercase).
    threshold : float
        Per-residue disorder score cutoff (default 0.5).

    Returns
    -------
    list[DisorderedRegion]
        List of contiguous disordered regions sorted by start position.
    """
    scores: list[float] = meta.predict_disorder(sequence)
    disordered_mask = [s >= threshold for s in scores]

    regions: list[DisorderedRegion] = []
    i = 0
    n = len(sequence)
    while i < n:
        if disordered_mask[i]:
            j = i
            while j < n and disordered_mask[j]:
                j += 1
            regions.append(DisorderedRegion(start=i, end=j, sequence=sequence[i:j]))
            i = j
        else:
            i += 1
    return regions


def find_lir_in_region(
    region: DisorderedRegion,
    motif_types: list[str] | None = None,
) -> Iterator[tuple[str, re.Match[str]]]:
    """
    Yield (motif_type, match) for every LIR hit inside *region*.

    Parameters
    ----------
    region : DisorderedRegion
    motif_types : list of str or None
        Subset of {"basic", "extended", "acidic_extended"}.
        Defaults to all three.
    """
    types = motif_types or list(_MOTIF_PATTERNS.keys())
    for mtype in types:
        pattern = _MOTIF_PATTERNS[mtype]
        for match in pattern.finditer(region.sequence):
            yield mtype, match


def analyze_protein(
    record: SeqRecord,
    threshold: float = 0.5,
    min_disorder_length: int = 1,
    motif_types: list[str] | None = None,
) -> ProteinResult:
    """
    Analyse one protein record for LIR motifs in disordered regions.

    Parameters
    ----------
    record : SeqRecord
    threshold : float
        Disorder score cutoff passed to ``predict_disorder``.
    min_disorder_length : int
        Minimum number of residues for a disordered region to be considered.
    motif_types : list[str] or None
        Which motif classes to search for (default: all).

    Returns
    -------
    ProteinResult
    """
    sequence = str(record._seq if hasattr(record, "_seq") else record.seq).upper()
    protein_id = record.id
    description = record.description

    # Run disorder prediction
    scores: list[float] = meta.predict_disorder(sequence)
    disordered_mask = [s >= threshold for s in scores]
    regions = [
        r for r in predict_disorder_from_mask(sequence, disordered_mask)
        if r.length >= min_disorder_length
    ]

    result = ProteinResult(
        protein_id=protein_id,
        protein_description=description,
        protein_length=len(sequence),
        n_disordered_regions=len(regions),
        n_disordered_residues=sum(r.length for r in regions),
    )

    for region in regions:
        region_scores = scores[region.start:region.end]
        for mtype, match in find_lir_in_region(region, motif_types):
            # match.group(1) because we used lookahead
            motif_seq = match.group(1)
            motif_local_start = match.start()
            motif_protein_start = region.start + motif_local_start  # 0-based
            motif_protein_end = motif_protein_start + len(motif_seq)  # exclusive

            # Mean disorder score over the motif residues
            motif_scores = scores[motif_protein_start:motif_protein_end]
            mean_score = sum(motif_scores) / len(motif_scores) if motif_scores else 0.0

            hit = MotifHit(
                protein_id=protein_id,
                protein_description=description,
                protein_length=len(sequence),
                motif_type=mtype,
                motif_sequence=motif_seq,
                motif_start_protein=motif_protein_start + 1,   # 1-based
                motif_end_protein=motif_protein_end,            # 1-based inclusive
                disorder_start_protein=region.start + 1,        # 1-based
                disorder_end_protein=region.end,                 # 1-based inclusive
                disorder_length=region.length,
                mean_disorder_score=round(mean_score, 4),
            )
            result.motif_hits.append(hit)

    return result


def predict_disorder_from_mask(
    sequence: str, mask: list[bool]
) -> list[DisorderedRegion]:
    """Build DisorderedRegion objects from a pre-computed boolean mask."""
    regions: list[DisorderedRegion] = []
    i, n = 0, len(sequence)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            regions.append(DisorderedRegion(start=i, end=j, sequence=sequence[i:j]))
            i = j
        else:
            i += 1
    return regions


def run_analysis(
    fasta_path: str,
    threshold: float = 0.5,
    min_disorder_length: int = 5,
    motif_types: list[str] | None = None,
) -> list[ProteinResult]:
    """
    Run the full pipeline on a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to input FASTA file.
    threshold : float
        Disorder score cutoff.
    min_disorder_length : int
        Minimum disordered region length to consider.
    motif_types : list[str] or None
        Motif classes to search for.

    Returns
    -------
    list[ProteinResult]
        One result object per protein in the FASTA.
    """
    results: list[ProteinResult] = []
    _iter = (
        SeqIO.parse(fasta_path, "fasta") if _BIOPYTHON else _parse_fasta(fasta_path)
    )
    for record in _iter:
        result = analyze_protein(
            record,
            threshold=threshold,
            min_disorder_length=min_disorder_length,
            motif_types=motif_types,
        )
        results.append(result)
    return results
