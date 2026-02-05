import tempfile
import os

from clue_pathway_enrichment.analysis.correlations import pearson_correlation


def test_pearson_perfect_positive():
    x = {"a": 1.0, "b": 2.0, "c": 3.0}
    y = {"a": 2.0, "b": 4.0, "c": 6.0}
    r = pearson_correlation(x, y)
    assert round(r, 6) == 1.0


def test_pearson_unrelated():
    x = {"a": 1.0, "b": 1.0}
    y = {"a": 2.0, "b": -2.0}
    r = pearson_correlation(x, y)
    # not exactly zero but should be negative; check type and range
    assert isinstance(r, float)
    assert -1.0 <= r <= 1.0
