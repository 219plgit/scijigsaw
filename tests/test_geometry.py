import numpy as np
import pytest

from scijigsaw import Piece, coverage, TAB, SOCKET
from scijigsaw.geometry import edge_path


def test_a_tab_protrudes_and_a_socket_indents():
    """The whole rule of the tool: complementary geometry."""
    flat = np.array(edge_path((0, 0), (1, 0)))
    tab = np.array(edge_path((0, 0), (1, 0), [(0.5, TAB)]))
    sock = np.array(edge_path((0, 0), (1, 0), [(0.5, SOCKET)]))
    # bottom edge travels left->right, so the outward normal points DOWN (-y)
    assert tab[:, 1].min() < flat[:, 1].min()    # tab pushes outward
    assert sock[:, 1].max() > flat[:, 1].max()   # socket cuts inward


def test_piece_polygon_is_closed_and_nondegenerate():
    p = Piece(0, 0, 2, 1).right(0.5, SOCKET).left(0.5, TAB)
    poly = p.polygon()
    assert poly.ndim == 2 and poly.shape[1] == 2
    assert len(poly) > 20
    assert np.isfinite(poly).all()


def test_ladder_places_only_engaged_subsites():
    """A partner brings a tab per contact it engages; the rest stay empty."""
    coords = [0.2, 0.4, 0.6, 0.8, 1.0]
    full = Piece(0, 0, 1, 2).ladder("right", coords, {0, 1, 2, 3, 4}, SOCKET)
    partial = Piece(0, 0, 1, 2).ladder("right", coords, {2, 3}, TAB)
    assert len(full.edges["right"]) == 5      # hub carries the whole site
    assert len(partial.edges["right"]) == 2   # partner engages 2 of 5


def test_coverage_is_a_ratio_not_an_affinity():
    assert coverage(5, 5) == 1.0
    assert coverage(2, 5) == 0.4
    with pytest.raises(ValueError):
        coverage(6, 5)


def test_degenerate_edge_rejected():
    with pytest.raises(ValueError, match="degenerate"):
        edge_path((0, 0), (0, 0))
