"""Shared pytest fixtures and mocks for LIR Motif Finder tests."""

from __future__ import annotations

from unittest.mock import patch
import pytest


# Mock metapredict so tests run without a real GPU / network download
@pytest.fixture(autouse=True)
def mock_metapredict():
    """Replace metapredict.predict_disorder with a deterministic stub.

    The stub returns a fixed disorder profile:
      - residues 5-14 are considered disordered (score 0.8)
      - all other residues are ordered (score 0.1)
    """
    def _fake_predict(sequence, *args, **kwargs):
        scores = []
        for i, _ in enumerate(sequence):
            scores.append(0.8 if 5 <= i < 15 else 0.1)
        return scores

    with patch("lirmotiffinder.core.meta.predict_disorder", side_effect=_fake_predict):
        yield
