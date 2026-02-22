"""
Unit tests for LIR Motif Finder.

Run with:
    pytest tests/ -v
"""

from __future__ import annotations

import re
import pytest

from lirmotiffinder.core import (
    DisorderedRegion,
    find_lir_in_region,
    predict_disorder_from_mask,
    _MOTIF_PATTERNS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def basic_lir_region():
    """A disordered region containing a basic LIR motif (WXXL)."""
    seq = "AADDEEWDDLKK"
    return DisorderedRegion(start=0, end=len(seq), sequence=seq)


@pytest.fixture
def extended_lir_region():
    """A disordered region containing an extended LIR motif (DWXXL)."""
    seq = "AADEWDDLKK"
    return DisorderedRegion(start=5, end=5 + len(seq), sequence=seq)


@pytest.fixture
def acidic_extended_region():
    """A disordered region containing an acidic-extended LIR motif (DDWXXL)."""
    seq = "AADDWDDIAA"
    return DisorderedRegion(start=10, end=10 + len(seq), sequence=seq)


@pytest.fixture
def no_motif_region():
    """A disordered region with no LIR motif."""
    seq = "AAKKKSSSNNN"
    return DisorderedRegion(start=0, end=len(seq), sequence=seq)


# ---------------------------------------------------------------------------
# Pattern tests
# ---------------------------------------------------------------------------
class TestMotifPatterns:
    def test_basic_matches_WxxL(self):
        assert _MOTIF_PATTERNS["basic"].search("WAAL")

    def test_basic_matches_FxxI(self):
        assert _MOTIF_PATTERNS["basic"].search("FAAI")

    def test_basic_matches_YxxV(self):
        assert _MOTIF_PATTERNS["basic"].search("YAAV")

    def test_basic_no_match_without_anchor(self):
        # First position must be W/F/Y
        assert not _MOTIF_PATTERNS["basic"].search("KDDL")

    def test_basic_no_match_without_Cterm(self):
        # Last position must be L/I/V
        assert not _MOTIF_PATTERNS["basic"].search("WAAK")

    def test_extended_requires_DE_prefix(self):
        assert _MOTIF_PATTERNS["extended"].search("DWAAL")
        assert _MOTIF_PATTERNS["extended"].search("EWAAL")
        # Without D/E prefix, should NOT match extended
        assert not _MOTIF_PATTERNS["extended"].search("KWAAL")

    def test_acidic_extended_requires_two_DE(self):
        assert _MOTIF_PATTERNS["acidic_extended"].search("DDWAAL")
        assert _MOTIF_PATTERNS["acidic_extended"].search("DEWAAL")
        assert _MOTIF_PATTERNS["acidic_extended"].search("EEWAAL")
        # Only one D/E should NOT match acidic_extended
        assert not _MOTIF_PATTERNS["acidic_extended"].search("DWAAL")

    def test_overlapping_motifs_found(self):
        """Lookahead allows overlapping matches."""
        seq = "WAALWDDL"
        matches = list(_MOTIF_PATTERNS["basic"].finditer(seq))
        assert len(matches) == 2


# ---------------------------------------------------------------------------
# DisorderedRegion tests
# ---------------------------------------------------------------------------
class TestDisorderedRegion:
    def test_length(self):
        r = DisorderedRegion(start=3, end=13, sequence="A" * 10)
        assert r.length == 10

    def test_sequence_stored(self):
        r = DisorderedRegion(start=0, end=5, sequence="ACDEF")
        assert r.sequence == "ACDEF"


# ---------------------------------------------------------------------------
# find_lir_in_region tests
# ---------------------------------------------------------------------------
class TestFindLirInRegion:
    def test_finds_basic_motif(self, basic_lir_region):
        hits = list(find_lir_in_region(basic_lir_region, motif_types=["basic"]))
        assert len(hits) >= 1
        motif_seqs = [m.group(1) for _, m in hits]
        assert any(re.match(r"[WFY][A-Z]{2}[LIV]", s) for s in motif_seqs)

    def test_finds_extended_motif(self, extended_lir_region):
        hits = list(find_lir_in_region(extended_lir_region, motif_types=["extended"]))
        assert len(hits) >= 1

    def test_finds_acidic_extended_motif(self, acidic_extended_region):
        hits = list(find_lir_in_region(acidic_extended_region, motif_types=["acidic_extended"]))
        assert len(hits) >= 1

    def test_no_hits_in_plain_region(self, no_motif_region):
        hits = list(find_lir_in_region(no_motif_region))
        assert hits == []

    def test_all_motif_types_by_default(self, basic_lir_region):
        """Default should search all motif types."""
        hits = list(find_lir_in_region(basic_lir_region))
        motif_types_found = {mtype for mtype, _ in hits}
        assert "basic" in motif_types_found

    def test_motif_type_filtering(self, basic_lir_region):
        hits_basic = list(find_lir_in_region(basic_lir_region, motif_types=["basic"]))
        hits_extended = list(find_lir_in_region(basic_lir_region, motif_types=["extended"]))
        # basic should have hits; extended may or may not
        all_types = {mtype for mtype, _ in hits_basic}
        assert all_types == {"basic"}


# ---------------------------------------------------------------------------
# predict_disorder_from_mask tests
# ---------------------------------------------------------------------------
class TestPredictDisorderFromMask:
    def test_single_region(self):
        seq = "ACDEFGHIJK"
        mask = [False] * 3 + [True] * 4 + [False] * 3
        regions = predict_disorder_from_mask(seq, mask)
        assert len(regions) == 1
        assert regions[0].start == 3
        assert regions[0].end == 7
        assert regions[0].sequence == seq[3:7]

    def test_multiple_regions(self):
        seq = "A" * 10
        mask = [True, True, False, False, True, True, True, False, False, False]
        regions = predict_disorder_from_mask(seq, mask)
        assert len(regions) == 2
        assert regions[0].length == 2
        assert regions[1].length == 3

    def test_all_ordered(self):
        seq = "ACDEF"
        mask = [False] * 5
        regions = predict_disorder_from_mask(seq, mask)
        assert regions == []

    def test_all_disordered(self):
        seq = "ACDEF"
        mask = [True] * 5
        regions = predict_disorder_from_mask(seq, mask)
        assert len(regions) == 1
        assert regions[0].length == 5
