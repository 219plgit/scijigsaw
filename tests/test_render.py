"""The renderer: an interaction table in, a jigsaw board out."""
import os

import pandas as pd
import pytest

from scijigsaw import Board, render
from scijigsaw.render import kd_to_coverage

HERE = os.path.dirname(__file__)
PROT = os.path.join(HERE, "..", "examples", "vamp2", "proteins.csv")
INTS = os.path.join(HERE, "..", "examples", "vamp2", "interactions.csv")


def test_kd_maps_to_contact_coverage():
    """Coverage proxies the EXTENT of interface engaged. It is not an affinity."""
    assert kd_to_coverage(0.5) == 5
    assert kd_to_coverage(50_000) == 2
    assert kd_to_coverage(None) is None


def test_render_writes_a_file(tmp_path):
    out = tmp_path / "board.png"
    b = Board(pd.read_csv(PROT), pd.read_csv(INTS))
    b.draw(str(out))
    assert out.exists() and out.stat().st_size > 10_000


def test_render_from_csv_end_to_end(tmp_path):
    out = tmp_path / "vamp2.pdf"
    b, path = render(PROT, INTS, str(out))
    assert os.path.exists(path)
    assert "VAMP2" in b.backbone


def test_report_names_the_bridges(tmp_path):
    b = Board(pd.read_csv(PROT), pd.read_csv(INTS))
    r = b.report()
    assert "Complexin" in r and "Syt-1" in r
    assert "bridges" in r


def test_a_partner_with_no_site_is_still_placed(tmp_path):
    """A missing site must not silently drop a protein from the board."""
    e = pd.read_csv(INTS)
    e.loc[e.protein_b == "Synaptophysin", ["site_on_a", "site_on_b"]] = None
    b = Board(pd.read_csv(PROT), e)
    placed = (set(b.backbone) | set(b.bridges) | set(b.pendants)
              | {c for c, _ in b.contenders})
    assert "Synaptophysin" in placed
