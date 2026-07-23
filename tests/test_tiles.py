"""The printable cut-out kit: keyed connectors and two variants.

These pin the connector semantics (a tab and a socket per interface, each with a
UNIQUE keyed shape so it fits only its partner), the student/teacher split (the
student set carries no solution labels), and that tiles do not collide on the
page -- a kit whose pieces overlap cannot be cut.
"""
import os

import pandas as pd
import pytest

from scijigsaw.geometry import KNOB, Piece, SOCKET, TAB
from scijigsaw.render import Board
from scijigsaw.tiles import TILE, A4_H, A4_W, TileKit, tiles

HERE = os.path.dirname(__file__)
VAMP2 = os.path.join(HERE, "..", "examples", "vamp2")


@pytest.fixture
def kit():
    b = Board(pd.read_csv(os.path.join(VAMP2, "proteins.csv")),
              pd.read_csv(os.path.join(VAMP2, "interactions.csv")))
    return TileKit(b)


def test_one_tile_per_protein(kit):
    assert set(kit.tiles) == set(kit.b.adj)


def test_every_interface_is_one_tab_and_one_socket(kit):
    assert len(kit.connectors) == len(kit.b.edges)
    for c in kit.connectors:
        assert c.tab != c.socket
        tab = [f for f in kit.tiles[c.tab].feats if f.code == c.code]
        soc = [f for f in kit.tiles[c.socket].feats if f.code == c.code]
        assert len(tab) == 1 and tab[0].kind == TAB
        assert len(soc) == 1 and soc[0].kind == SOCKET
        assert tab[0].cov == soc[0].cov == c.cov


def test_connectors_are_shape_keyed(kit):
    """A tab and its socket share a key; different interfaces differ, so a tab
    fits only its true partner."""
    for c in kit.connectors:
        tab = next(f for f in kit.tiles[c.tab].feats if f.code == c.code)
        soc = next(f for f in kit.tiles[c.socket].feats if f.code == c.code)
        assert (tab.shape, tab.size) == (soc.shape, soc.size) == (c.shape, c.size)
    keys = [(c.shape, c.size) for c in kit.connectors]
    assert len(set(keys)) == len(keys), "connector silhouettes must be distinct"


def test_bridge_and_alt_are_marked(kit):
    kinds = {c.kind for c in kit.connectors}
    assert "bridge" in kinds and "alt" in kinds
    bridges = {(c.tab, c.socket) for c in kit.connectors if c.kind == "bridge"}
    assert ("VAMP2", "Complexin") in bridges
    assert ("Syntaxin-1A", "Complexin") in bridges
    alts = [c for c in kit.connectors if c.kind == "alt"]
    assert alts and all("either/or" in c.note for c in alts)


def test_tiles_do_not_overlap_on_a_page(kit):
    positions, per = kit._slots()
    ordered = [kit.tiles[n] for n in kit._pieces() if n in kit.tiles][:per]

    def bbox(t, cx, cy):
        pc = Piece(cx, cy, TILE, TILE)
        for f in t.feats:
            coord = (cy + f.frac * TILE) if f.edge in ("right", "left") \
                else (cx + f.frac * TILE)
            getattr(pc, f.edge)(coord, f.kind, f.size, f.shape)
        p = pc.polygon()
        return (p[:, 0].min() - 8, p[:, 0].max() + 8,
                p[:, 1].min() - 9, p[:, 1].max() + 9)

    boxes = [bbox(t, x, y) for (x, y), t in zip(positions, ordered)]
    for x0, x1, y0, y1 in boxes:
        assert x0 > 0 and x1 < A4_W and y0 > 0 and y1 < A4_H
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            disjoint = a[1] < b[0] or b[1] < a[0] or a[3] < b[2] or b[3] < a[2]
            assert disjoint, "tiles overlap; the kit would not cut cleanly"


def test_renders_both_variants(tmp_path, kit):
    for v in ("student", "teacher"):
        for ext in ("pdf", "svg"):
            out = str(tmp_path / f"kit_{v}.{ext}")
            kit.draw(out, variant=v)
            assert os.path.exists(out) and os.path.getsize(out) > 2000


def test_student_hides_the_solution(tmp_path, kit):
    stu = str(tmp_path / "s.svg"); tea = str(tmp_path / "t.svg")
    kit.draw(stu, variant="student"); kit.draw(tea, variant="teacher")
    s, t = open(stu).read(), open(tea).read()
    assert "VAMP2" in s and "VAMP2" in t          # names on both
    for label in ("backbone", "seat after", "BRIDGE"):
        assert label not in s                      # no solution on the class set
        assert label in t                          # present on the answer key


def test_public_entry_point(tmp_path):
    _, both = tiles(os.path.join(VAMP2, "proteins.csv"),
                    os.path.join(VAMP2, "interactions.csv"),
                    str(tmp_path / "kit.pdf"), variant="both")
    assert len(both) == 2 and all(os.path.exists(p) for p in both)
    _, one = tiles(os.path.join(VAMP2, "proteins.csv"),
                   os.path.join(VAMP2, "interactions.csv"),
                   str(tmp_path / "solo.pdf"), variant="teacher")
    assert len(one) == 1 and os.path.exists(one[0])
