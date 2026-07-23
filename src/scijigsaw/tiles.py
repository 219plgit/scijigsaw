"""Printable cut-out tile kit: a board you can build with scissors.

`scijigsaw-render` draws the board already assembled. This module does the
opposite. It lays every protein out as a SEPARATE tile on A4 pages, so the sheet
can be printed, cut along the solid outlines, and assembled by hand.

Every interface is drawn as a complementary pair -- a *tab* (protrudes) and a
*socket* (indents) -- and, crucially, EACH interface is given its own connector
SHAPE (a keyed silhouette: round / wedge / square, in four sizes). Complementary
pieces therefore fit only their true partner, which is what makes the kit
self-correcting: a piece that does not fit signals an encoding it does not match
-- Montessori's *control of error*, in cardboard.

Two variants are produced from the same geometry:

  STUDENT   Each tile carries ONLY the protein name. Learners build the complex
            by matching connector shapes and reasoning about the biology.
  TEACHER   The answer key. Adds the connector number on every tab/socket, the
            interface contact coverage n/N, and the precedence (bridge, "seat
            last") and exclusion (dashed, "either/or") cues, plus a key page.

Because both variants share identical cut geometry, a teacher's solved model and
a student's tiles are the same physical pieces. Topology (backbone / bridge /
pendant / alternative occupancy), colour, and evolutionary rings are reused from
:class:`render.Board`.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Polygon, Rectangle, RegularPolygon

from .geometry import SOCKET, TAB, Piece
from .render import AGE_RINGS, INK, MUTED, PALETTE, STOP, Board

# ------------------------------------------------------------------ page metrics
A4_W, A4_H = 210.0, 297.0          # mm
MARGIN = 12.0
TILE = 46.0                        # side of a tile, mm
CLEAR = 16.0                       # clearance for tabs + badges around a tile
BADGE_R = 4.4
N_SUB = 5

# connector keying: a unique (shape, size) per interface, so a tab fits only its
# own socket. Twelve distinct keys before any repeat.
SHAPES = ["round", "wedge", "square"]
SIZES = [0.085, 0.110, 0.135, 0.160]


def _key(code: int):
    i = code - 1
    return SHAPES[i % 3], SIZES[(i // 3) % len(SIZES)]


@dataclass
class Connector:
    code: int
    tab: str
    socket: str
    cov: int
    kind: str
    note: str = ""
    shape: str = "round"
    size: float = 0.11


@dataclass
class _Feat:
    code: int
    kind: int          # TAB or SOCKET (geometry sign)
    ckind: str         # connector kind
    partner: str
    cov: int
    shape: str = "round"
    size: float = 0.11
    edge: str = ""
    frac: float = 0.5


@dataclass
class _Tile:
    name: str
    role: str
    feats: List[_Feat] = field(default_factory=list)
    banner: str = ""


class TileKit:
    """Turns a :class:`render.Board` into a printable set of cut-out tiles."""

    def __init__(self, board: Board):
        self.b = board
        self.connectors: List[Connector] = []
        self.tiles: Dict[str, _Tile] = {}
        self._build()

    # ------------------------------------------------------------- assignment
    def _role(self, name: str) -> str:
        if name in self.b.bridges:
            return "bridge"
        if any(c == name for c, _ in self.b.contenders):
            return "alt"
        if name in self.b.pendants:
            return "pendant"
        return "backbone"

    def _classify(self, a: str, b: str):
        for br, parents in self.b.bridges.items():
            if br in (a, b) and (a in parents or b in parents):
                parent = b if a == br else a
                if parent in parents:
                    return "bridge", parent, br, ""
        for c, (host, site, winner) in self.b.contenders:
            if {a, b} == {c, host}:
                return "alt", c, host, f"either/or with {winner}"
        if self.b.pendants.get(a) == b:
            return "pendant", a, b, ""
        if self.b.pendants.get(b) == a:
            return "pendant", b, a, ""
        order = self.b.backbone

        def rank(n):
            return order.index(n) if n in order else len(order) + sorted(self.b.adj).index(n)
        socket_owner, tab_owner = (a, b) if rank(a) <= rank(b) else (b, a)
        return "backbone", tab_owner, socket_owner, ""

    def _build(self):
        for n in self._pieces():
            self.tiles[n] = _Tile(n, self._role(n))

        code = 0
        seen = set()
        for r in self.b.edges.itertuples():
            a, b = r.protein_a, r.protein_b
            key = frozenset((a, b))
            if key in seen or a not in self.tiles or b not in self.tiles:
                continue
            seen.add(key)
            kind, tab_owner, socket_owner, note = self._classify(a, b)
            cov = self.b._cov(a, b)
            code += 1
            shape, size = _key(code)
            self.connectors.append(
                Connector(code, tab_owner, socket_owner, cov, kind, note, shape, size))
            self.tiles[tab_owner].feats.append(
                _Feat(code, TAB, kind, socket_owner, cov, shape, size))
            self.tiles[socket_owner].feats.append(
                _Feat(code, SOCKET, kind, tab_owner, cov, shape, size))

        for name, t in self.tiles.items():
            if t.role == "bridge":
                p, q = self.b.bridges[name]
                t.banner = f"BRIDGE \u00b7 seat after {p} + {q}"
            elif t.role == "alt":
                for c, (host, site, winner) in self.b.contenders:
                    if c == name:
                        t.banner = ("on %s\u2019s %s\n\u2298 excludes %s"
                                    % (host, site.replace("_", " "), winner))
            self._place(t)

    def _pieces(self) -> List[str]:
        order = list(self.b.backbone)
        order += sorted(self.b.bridges)
        order += sorted(self.b.pendants)
        order += sorted(c for c, _ in self.b.contenders)
        for n in sorted(self.b.adj):
            if n not in order:
                order.append(n)
        return order

    @staticmethod
    def _place(t: _Tile):
        edges = ["right", "left", "top", "bottom"]
        buckets: Dict[str, List[_Feat]] = defaultdict(list)
        for i, f in enumerate(t.feats):
            buckets[edges[i % 4]].append(f)
        for e, items in buckets.items():
            m = len(items)
            for j, f in enumerate(items):
                f.edge = e
                f.frac = 0.24 + 0.52 * ((j + 1) / (m + 1))

    # ----------------------------------------------------------------- drawing
    def _style(self, n):
        return self.b._style(n)

    def _age(self, n):
        return self.b._age(n)

    @staticmethod
    def _rings(ax, x, y, n, col=INK, r0=1.6):
        ax.add_patch(Circle((x, y), r0 * (n + 1), facecolor="white",
                            edgecolor=col, lw=0.6, zorder=9))
        for k in range(1, n + 1):
            ax.add_patch(Circle((x, y), r0 * k, facecolor="none",
                                edgecolor=col, lw=0.7, zorder=10))

    def _draw_tile(self, ax, t: _Tile, cx: float, cy: float, variant: str):
        S = TILE
        pc = Piece(cx, cy, S, S)
        badges = []
        for f in t.feats:
            if f.edge == "right":
                Y = cy + f.frac * S
                pc.right(Y, f.kind, f.size, f.shape)
                badges.append((f, cx + S + 6.5, Y))
            elif f.edge == "left":
                Y = cy + f.frac * S
                pc.left(Y, f.kind, f.size, f.shape)
                badges.append((f, cx - 6.5, Y))
            elif f.edge == "top":
                X = cx + f.frac * S
                pc.top(X, f.kind, f.size, f.shape)
                badges.append((f, X, cy + S + 6.5))
            else:
                X = cx + f.frac * S
                pc.bottom(X, f.kind, f.size, f.shape)
                badges.append((f, X, cy - 6.5))

        fc, ec = self._style(t.name)
        dashed = t.role == "alt" and variant == "teacher"
        ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                             edgecolor=INK, lw=1.7,
                             linestyle=(0, (4, 2)) if dashed else "solid",
                             joinstyle="round", zorder=3))

        fs = 9.5 if len(t.name) <= 9 else (8.0 if len(t.name) <= 13 else 6.8)
        yname = cy + S / 2 + (1.5 if variant == "teacher" else 0.0)
        ax.text(cx + S / 2, yname, t.name, color=INK, fontsize=fs,
                ha="center", va="center", fontweight="bold", zorder=6)
        self._rings(ax, cx + 5.5, cy + S - 5.5, self._age(t.name))

        if variant == "teacher":
            ax.text(cx + S / 2, cy + S / 2 - 4.6, f"[{t.role}]", color=ec,
                    fontsize=5.6, ha="center", va="center", zorder=6)
            if t.banner:
                col = STOP if t.role == "alt" else "#3d7ab8"
                ax.text(cx + S / 2, cy + 4.4, t.banner, color=col, fontsize=5.2,
                        ha="center", va="center", fontweight="bold", zorder=6,
                        linespacing=0.95)
            for f, bx, by in badges:
                self._badge(ax, f, bx, by, ec)

    def _badge(self, ax, f: _Feat, x, y, ec):
        is_tab = f.kind == TAB
        if f.ckind == "bridge":
            ax.add_patch(RegularPolygon((x, y), 6, radius=BADGE_R + 0.6,
                         orientation=0.52,
                         facecolor=ec if is_tab else "white",
                         edgecolor=ec, lw=1.0, zorder=11))
        else:
            ax.add_patch(Circle((x, y), BADGE_R,
                         facecolor=ec if is_tab else "white",
                         edgecolor=ec, lw=1.0, zorder=11))
        ax.text(x, y, str(f.code), color="white" if is_tab else ec,
                fontsize=5.6, ha="center", va="center", fontweight="bold", zorder=12)
        arrow = "\u25b6" if is_tab else "\u25c0"
        mark = " \u2205" if f.ckind == "alt" else (" \u2229" if f.ckind == "bridge" else "")
        _nn = f"{f.cov}/{N_SUB}" if f.cov is not None else "\u2014"
        ax.text(x, y - BADGE_R - 2.0, f"{arrow}{_nn}{mark}", color=ec,
                fontsize=4.6, ha="center", va="center", zorder=12)

    # --------------------------------------------------------------- key pages
    def _draw_key(self, ax):
        ax.text(MARGIN, A4_H - MARGIN - 2, "Assembly key (teacher)", color=INK,
                fontsize=13, fontweight="bold", va="top")
        steps = [
            "1.  Print at 100% (no 'fit to page'), then cut along the solid outlines.",
            "2.  Every interface has its own connector shape: a tab fits ONLY its matching socket.",
            "3.  On this teacher set the number confirms it: tab \u25b6 n meets socket \u25c0 n.",
            "4.  Build the backbone first: " + " - ".join(self.b.backbone) + ".",
            "5.  Seat a bridge \u2229 only after BOTH partners beneath it are placed (precedence).",
            "6.  A dashed tile \u2205 shares its rival's socket: place it OR its rival, never both.",
        ]
        y = A4_H - MARGIN - 12
        for s in steps:
            ax.text(MARGIN, y, s, color=INK, fontsize=7.4, va="top")
            y -= 6.6
        y -= 3
        ax.text(MARGIN, y, "Connectors", color=INK, fontsize=9, fontweight="bold", va="top")
        y -= 7
        hdr = f"{'#':>2}   {'tab \u25b6':<16}{'socket \u25c0':<16}{'n/N':<6}{'shape':<8}{'type':<10}note"
        ax.text(MARGIN, y, hdr, color=MUTED, fontsize=6.6, va="top", family="monospace")
        y -= 5
        for c in self.connectors:
            line = (f"{c.code:>2}   {c.tab:<16}{c.socket:<16}{c.cov}/{N_SUB}   "
                    f"{c.shape:<8}{c.kind:<10}{c.note}")
            ax.text(MARGIN, y, line, color=INK, fontsize=6.6, va="top", family="monospace")
            y -= 5
            if y < MARGIN + 30:
                break
        ax.text(MARGIN, MARGIN + 16,
                "Counts of the assembly orders permitted by this board are printed by "
                "`scijigsaw-count`. The kit encodes the same precedence (bridges) and exclusion "
                "(dashed tiles) that the enumerator counts; it does not assert a biological pathway.",
                color=MUTED, fontsize=6.2, va="top")
        ax.text(MARGIN, MARGIN + 4,
                "scijigsaw \u00b7 interface geometry as a constraint on protein-assembly order",
                color=MUTED, fontsize=6.0, va="top", style="italic")

    def _draw_student_intro(self, ax):
        ax.text(MARGIN, A4_H - MARGIN - 2, "Build the complex", color=INK,
                fontsize=13, fontweight="bold", va="top")
        steps = [
            "1.  Print at 100% (no 'fit to page'), then cut out every tile along its outline.",
            "2.  Each tile is one protein. Its bumps (tabs) and notches (sockets) are its encoded interfaces.",
            "3.  A tab fits ONLY the socket of the SAME shape and size \u2014 that is its true partner.",
            "4.  A piece with two bumps is a bridge: it sits across a seam and can go in only once",
            "     both proteins beneath it are already joined.",
            "5.  If two pieces want the same notch, only one can take it: they are alternatives.",
            "6.  When every tab has found its socket and nothing is forced, you have built the complex.",
        ]
        y = A4_H - MARGIN - 12
        for s in steps:
            ax.text(MARGIN, y, s, color=INK, fontsize=8.2, va="top")
            y -= 8.0
        ax.text(MARGIN, MARGIN + 6,
                "scijigsaw \u00b7 a jigsaw where the pieces refuse to fit unless the interfaces agree",
                color=MUTED, fontsize=6.4, va="top", style="italic")

    # ------------------------------------------------------------------ layout
    def _slots(self):
        usable_w = A4_W - 2 * MARGIN
        cols = max(1, min(2, int(usable_w // (TILE + 2 * CLEAR))))
        step_x = usable_w / cols
        top = A4_H - MARGIN - 20
        step_y = TILE + 2 * CLEAR
        rows = max(1, int((top - MARGIN) // step_y))
        positions = []
        for r in range(rows):
            for c in range(cols):
                x = MARGIN + c * step_x + (step_x - TILE) / 2
                y = top - (r + 1) * step_y + (step_y - TILE) / 2
                positions.append((x, y))
        return positions, rows * cols

    def draw(self, out: str, variant: str = "teacher", title: str = None):
        if variant not in ("teacher", "student"):
            raise ValueError("variant must be 'teacher' or 'student'")
        title = title or ("Scientific Jigsaw \u2014 cut-out kit "
                           + ("(answer key)" if variant == "teacher" else "(class set)"))
        tiles = [self.tiles[n] for n in self._pieces() if n in self.tiles]
        positions, per_page = self._slots()

        def header(ax, page, npages):
            ax.text(MARGIN, A4_H - MARGIN - 2, title, color=INK, fontsize=12.5,
                    fontweight="bold", va="top")
            ax.text(A4_W - MARGIN, A4_H - MARGIN - 2, f"sheet {page}/{npages}",
                    color=MUTED, fontsize=7.5, ha="right", va="top")
            if variant == "teacher":
                sub = ("Cut on solid lines. Each interface is a unique connector shape; the number "
                       "confirms the match. \u2229 bridge = seat last. \u2205 dashed = either/or.")
            else:
                sub = ("Cut on solid lines. A bump (tab) fits only the notch (socket) of the same "
                       "shape and size \u2014 its true partner.")
            ax.text(MARGIN, A4_H - MARGIN - 10, sub, color=MUTED, fontsize=6.8, va="top")

        n_tile_pages = (len(tiles) + per_page - 1) // per_page
        npages = n_tile_pages + 1

        def render_page(pdf, pg):
            fig = plt.figure(figsize=(A4_W / 25.4, A4_H / 25.4), facecolor="white")
            ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, A4_W); ax.set_ylim(0, A4_H)
            ax.set_aspect("equal"); ax.axis("off")
            header(ax, pg + 1, npages)
            for (x, y), t in zip(positions, tiles[pg * per_page:(pg + 1) * per_page]):
                self._draw_tile(ax, t, x, y, variant)
            pdf.savefig(fig, facecolor="white"); plt.close(fig)

        if out.lower().endswith(".pdf"):
            with PdfPages(out) as pdf:
                for pg in range(n_tile_pages):
                    render_page(pdf, pg)
                fig = plt.figure(figsize=(A4_W / 25.4, A4_H / 25.4), facecolor="white")
                ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, A4_W); ax.set_ylim(0, A4_H)
                ax.set_aspect("equal"); ax.axis("off")
                if variant == "teacher":
                    self._draw_key(ax)
                else:
                    render_legend(ax)
                pdf.savefig(fig, facecolor="white"); plt.close(fig)
        else:
            cols = 2
            rows_needed = (len(tiles) + cols - 1) // cols
            H = MARGIN * 2 + 20 + rows_needed * (TILE + 2 * CLEAR)
            fig = plt.figure(figsize=(A4_W / 25.4, H / 25.4), facecolor="white")
            ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, A4_W); ax.set_ylim(0, H)
            ax.set_aspect("equal"); ax.axis("off")
            ax.text(MARGIN, H - MARGIN - 2, title, color=INK, fontsize=12.5,
                    fontweight="bold", va="top")
            step_x = (A4_W - 2 * MARGIN) / cols
            step_y = TILE + 2 * CLEAR
            top = H - MARGIN - 16
            for i, t in enumerate(tiles):
                c, r = i % cols, i // cols
                x = MARGIN + c * step_x + (step_x - TILE) / 2
                y = top - (r + 1) * step_y + (step_y - TILE) / 2
                self._draw_tile(ax, t, x, y, variant)
            fig.savefig(out, facecolor="white", bbox_inches="tight"); plt.close(fig)
        return out


def tiles(proteins_csv, interactions_csv, out, variant="both", title=None):
    """Build printable cut-out kit(s) from an interaction table.

    variant: 'teacher' (answer key), 'student' (names only) or 'both'. For
    'both', '_teacher'/'_student' are inserted before the extension.
    """
    board = Board(pd.read_csv(proteins_csv), pd.read_csv(interactions_csv))
    kit = TileKit(board)
    outs = []
    wanted = ["teacher", "student"] if variant == "both" else [variant]
    for v in wanted:
        if variant == "both":
            stem, _, ext = out.rpartition(".")
            path = f"{stem}_{v}.{ext}" if ext else f"{out}_{v}"
        else:
            path = out
        kit.draw(path, variant=v, title=title)
        outs.append(path)
    return kit, outs


_FUNC_LABELS = [
    ("fusion", "fusion / SNARE core"),
    ("regulation", "regulation"),
    ("retrieval", "endocytic retrieval"),
    ("chaperone", "chaperone / SM protein"),
    ("contact", "membrane contact"),
    ("other", "other / unspecified"),
]
_AGE_LABELS = [(3, "ancient — conserved since early eukaryotes (yeast homologue)"),
               (2, "metazoan — conserved across animals"),
               (1, "vertebrate-only — recent, least conserved")]


def _mini(ax, x, y, kind, shape, size, s=13.0, fc="#d7dbe2", ec=INK, dashed=False):
    """A small sample piece with one connector on its top edge."""
    pc = Piece(x, y, s, s)
    pc.top(x + s * 0.5, kind, size, shape)
    ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec,
                         lw=1.4, joinstyle="round",
                         linestyle=(0, (3, 2)) if dashed else "solid", zorder=3))


def render_legend(ax):
    """Draw the reference legend as a SINGLE COLUMN: sample on the left, text on
    the right, one row at a time down the page, so nothing can overlap."""
    x = MARGIN
    tx = x + 26.0            # text starts to the right of every sample
    y = A4_H - MARGIN - 4

    ax.text(x, y, "How to read the pieces", color=INK, fontsize=16,
            fontweight="bold", va="top")
    y -= 9
    ax.text(x, y, "This page decodes the four visual channels and two seating "
            "constraints.", color=MUTED,
            fontsize=9, va="top")
    y -= 13

    def section(title):
        nonlocal y
        ax.text(x, y, title, color=INK, fontsize=11, fontweight="bold", va="top")
        y -= 10

    def row(draw, text, h=13.0, size=8.5):
        nonlocal y
        cy = y - h / 2.0
        if draw is not None:
            draw(cy)
        ax.text(tx, cy, text, color=INK, fontsize=size, va="center")
        y -= h

    section("Shape — the encoded interfaces")
    row(lambda cy: _mini(ax, x, cy - 6.5, TAB, "round", 0.16),
        "Tab (bump): an interface this protein OFFERS.", h=16)
    row(lambda cy: _mini(ax, x, cy - 6.5, SOCKET, "round", 0.16),
        "Socket (notch): an interface it RECEIVES.", h=16)

    def _shapes(cy):
        for i, sh in enumerate(SHAPES):
            _mini(ax, x + i * 8.5, cy - 4.0, TAB, sh, 0.16, s=8.0)
    row(_shapes,
        "Each interface has a symbolic shape and size; a tab fits only its "
        "matching socket (a key, not a molecular surface).", h=15)
    y -= 2

    section("Colour — functional class")
    for fkey, lab in _FUNC_LABELS:
        fc, ec = PALETTE[fkey]
        row(lambda cy, fc=fc, ec=ec: ax.add_patch(
                Rectangle((x, cy - 3), 8, 6, facecolor=fc, edgecolor=ec, lw=1.0)),
            lab, h=9.4, size=8.0)
    y -= 2

    section("Rings — conservation tier")
    ax.text(x, y, "More rings = a more deeply conserved phylogenetic tier "
            "(annotation only; it does not affect placement or counting).",
            color=MUTED, fontsize=8, va="top")
    y -= 8
    for n, lab in _AGE_LABELS:
        row(lambda cy, n=n: TileKit._rings(ax, x + 4, cy, n, r0=1.5), lab,
            h=12, size=8.0)
    y -= 2

    section("On the teacher set only — the answer")

    def _numbadge(cy):
        ax.add_patch(Circle((x + 4, cy), BADGE_R, facecolor=PALETTE["fusion"][1],
                     edgecolor=PALETTE["fusion"][1], lw=1.0))
        ax.text(x + 4, cy, "3", color="white", fontsize=7, ha="center",
                va="center", fontweight="bold")
    row(_numbadge, "Connector number: tab 3 and socket 3 are the same interface "
        "— partners.", h=12)
    row(lambda cy: ax.text(x + 4, cy, "▶ n/N", color=INK, fontsize=8.6,
        ha="center", va="center"),
        "n/N (structure-derived only): hub residues contacted / overlap-group "
        "footprint; blank for curated relations. Not an order-count input.", h=12)
    row(lambda cy: ax.add_patch(RegularPolygon((x + 4, cy), 6, radius=BADGE_R + 0.6,
        orientation=0.52, facecolor="white", edgecolor=INK, lw=1.0)),
        "Bridge (∩): a two-socket piece. Seat it only after BOTH partners "
        "beneath are down (precedence).", h=12)
    row(lambda cy: _mini(ax, x, cy - 5.5, SOCKET, "round", 0.16, s=10,
        fc=PALETTE["retrieval"][0], dashed=True),
        "Dashed outline (∅): alternative occupancy. Place this piece OR its "
        "rival, never both (exclusion).", h=15)

    ax.text(x, MARGIN, "scijigsaw · interface geometry as a constraint on "
            "protein-assembly order", color=MUTED, fontsize=7.5, va="bottom",
            style="italic")

def legend(out, title=None):
    """Render the standalone legend page to a vector file."""
    fig = plt.figure(figsize=(A4_W / 25.4, A4_H / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, A4_W); ax.set_ylim(0, A4_H)
    ax.set_aspect("equal"); ax.axis("off")
    render_legend(ax)
    fig.savefig(out, facecolor="white"); plt.close(fig)
    return out
