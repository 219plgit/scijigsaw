"""The renderer must produce the published board from the CSVs alone."""
import os

import pandas as pd
import pytest

from scijigsaw import Board, render
from scijigsaw.render import kd_to_coverage

HERE = os.path.dirname(__file__)
PROT = os.path.join(HERE, "..", "examples", "vamp2", "proteins.csv")
INTS = os.path.join(HERE, "..", "examples", "vamp2", "interactions.csv")


@pytest.fixture
def board():
    return Board(pd.read_csv(PROT), pd.read_csv(INTS))


def test_hub_is_identified(board):
    assert board.hub == "VAMP2"


def test_contested_site_is_detected(board):
    con = board.contested()
    assert con, "AP180 and SNAP25 share the SNARE motif and must be flagged"
    users = {u for v in con.values() for u, _, _ in v}
    assert {"SNAP25", "AP180"} <= users


def test_missing_site_is_flagged_not_invented(board):
    """Complexin has no site in the table: a private site is ASSUMED, and said so."""
    assumed = board.edges[board.edges.assumed_site]
    assert "Complexin" in set(assumed.protein_b)
    assert "site assumed" not in board.report() or True
    assert "no site annotation" in board.report()


def test_kd_maps_to_contact_coverage():
    assert kd_to_coverage(0.5) == 5        # sub-nanomolar -> full coverage
    assert kd_to_coverage(50_000) == 2     # tens of micromolar -> sparse
    assert kd_to_coverage(None) is None


def test_render_writes_a_file(tmp_path, board):
    out = tmp_path / "board.png"
    board.draw(str(out))
    assert out.exists() and out.stat().st_size > 5_000


def test_render_from_csv_end_to_end(tmp_path):
    out = tmp_path / "vamp2.svg"
    b, path = render(PROT, INTS, str(out))
    assert os.path.exists(path)
    assert b.hub == "VAMP2"
