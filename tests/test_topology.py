"""The board's TOPOLOGY must be derived from the table, not hard-coded.

The paper claims Figure 2 is produced by the renderer from the deposited CSVs. These
tests are what make that claim checkable.
"""
import os

import pandas as pd
import pytest

from scijigsaw import Board

HERE = os.path.dirname(__file__)
EX = os.path.join(HERE, "..", "examples", "vamp2")


@pytest.fixture
def board():
    return Board(pd.read_csv(os.path.join(EX, "proteins.csv")),
                 pd.read_csv(os.path.join(EX, "interactions.csv")))


def test_backbone_is_the_snare_chain(board):
    assert set(board.backbone) == {"VAMP2", "SNAP25", "Syntaxin-1A"}
    # SNAP25 must lie between VAMP2 and syntaxin
    assert board.backbone.index("SNAP25") == 1


def test_bridges_are_derived_not_declared(board):
    """A bridge binds two ADJACENT proteins, so it cannot be seated until the seam
    between them closes. This is the constraint that does the combinatorial work."""
    assert set(board.bridges) == {"Complexin", "Syt-1"}
    assert set(board.bridges["Complexin"]) == {"VAMP2", "SNAP25"}
    assert set(board.bridges["Syt-1"]) == {"SNAP25", "Syntaxin-1A"}


def test_a_shared_site_is_not_always_competition(board):
    """Complexin's central helix contacts BOTH VAMP2 and SNAP25 -- one interface, two
    partners. They are adjacent, so complexin lies between them: a bridge, not a
    contested site. Getting this wrong benches half the board."""
    benched = {c for c, _ in board.contenders}
    assert "Complexin" not in benched
    assert "SNAP25" not in benched


def test_competition_requires_non_adjacent_partners(board):
    """AP180 and SNAP25 both take VAMP2's SNARE motif, and AP180 does not bind
    SNAP25 -- so one excludes the other."""
    benched = {c for c, _ in board.contenders}
    assert benched == {"AP180", "CALM"}
    for c, (host, site, winner) in board.contenders:
        assert host == "VAMP2" and site == "SNARE_motif" and winner == "SNAP25"


def test_pendants_bind_exactly_one_backbone_protein(board):
    assert board.pendants == {"Munc18-1": "Syntaxin-1A",
                              "SNCA": "VAMP2",
                              "Synaptophysin": "VAMP2"}


def test_every_protein_is_placed(board):
    placed = (set(board.backbone) | set(board.bridges) | set(board.pendants)
              | {c for c, _ in board.contenders})
    assert placed == set(board.adj), "a protein in the table must appear on the board"


def test_figure2_is_reproducible(tmp_path, board):
    out = tmp_path / "Figure2.pdf"
    board.draw(str(out))
    assert out.exists() and out.stat().st_size > 10_000
