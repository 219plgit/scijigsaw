"""Board renderer: an interaction table in, a jigsaw board out.

    scijigsaw-render examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
        --out board.pdf

This renders a GENERAL interaction graph, not merely a hub with spokes. The paper's
central mechanic is the BRIDGE PIECE -- a partner binding two proteins that are
themselves adjacent, and which therefore cannot be seated until the seam between them
is closed. A renderer that cannot draw a bridge cannot express the argument, so the
topology is DERIVED from the table rather than hard-coded:

    BACKBONE   the longest chain of mutually binding proteins. Drawn interlocking.
    BRIDGE     a piece binding two ADJACENT backbone proteins. Drawn above, spanning
               the seam, carrying a socket over each. This is precisely the precedence
               constraint that assembly.py counts.
    PENDANT    a piece binding exactly one backbone protein.
    CONTENDER  a piece whose site on a partner is already occupied. Benched -- overlap
               establishes NON-SIMULTANEITY, not competition; an enzyme may drive
               succession instead.

Encoding: shape = what binds where; colour = function; rings = evolutionary age;
n/N = interface contact coverage (extent of interface engaged, NOT affinity).
"""
from __future__ import annotations

import itertools
import math
from collections import defaultdict
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Polygon

from .geometry import KNOB, SOCKET, TAB, Piece

INK, MUTED, STOP = "#1a1d24", "#5b6675", "#c0392b"
PALETTE: Dict[str, tuple] = {
    "fusion":     ("#d7dbe2", "#78828f"),
    "regulation": ("#cfe0f2", "#3d7ab8"),
    "retrieval":  ("#cbe9e3", "#1f7a6f"),
    "chaperone":  ("#fbdc8a", "#a8760a"),
    "contact":    ("#ded3f0", "#5b3f8c"),
    "other":      ("#e6e8ec", "#8b95a7"),
}
AGE_RINGS = {"ancient": 3, "metazoan": 2, "vertebrate": 1}
N_SUB = 5
KD_LADDER = [(1.0, 5), (100.0, 4), (1_000.0, 3), (100_000.0, 2), (math.inf, 1)]


def kd_to_coverage(kd):
    """Kd -> subsites engaged. A DECLARED CONVENTION: contact coverage proxies the
    extent of interface engaged, not affinity."""
    if kd is None or (isinstance(kd, float) and math.isnan(kd)):
        return None
    for limit, n in KD_LADDER:
        if float(kd) < limit:
            return n
    return 1


def _rings(ax, x, y, n, col=INK, r0=0.042, z=9):
    ax.add_patch(plt.Circle((x, y), r0 * (n + 1), facecolor="white", edgecolor=col,
                            lw=0.7, zorder=z))
    for k in range(1, n + 1):
        ax.add_patch(plt.Circle((x, y), r0 * k, facecolor="none", edgecolor=col,
                                lw=0.8, zorder=z + 1))


class Board:
    def __init__(self, proteins: pd.DataFrame, interactions: pd.DataFrame):
        self.meta = proteins.set_index("name")
        e = interactions.copy()
        e["coverage"] = e["kd_nM"].apply(kd_to_coverage) if "kd_nM" in e else None
        self.edges = e

        self.adj: Dict[str, set] = defaultdict(set)
        self.cov: Dict[frozenset, int] = {}
        self.site: Dict[tuple, str] = {}
        for r in e.itertuples():
            self.adj[r.protein_a].add(r.protein_b)
            self.adj[r.protein_b].add(r.protein_a)
            self.cov[frozenset((r.protein_a, r.protein_b))] = r.coverage
            self.site[(r.protein_a, r.protein_b)] = getattr(r, "site_on_a", None)
            self.site[(r.protein_b, r.protein_a)] = getattr(r, "site_on_b", None)

        # The CORE comes first. A four-helix bundle such as the SNARE complex is a
        # CLIQUE, not a chain -- every helix contacts every other -- so a rule of the
        # form "its partners are mutually adjacent, therefore it bridges them" fires on
        # every member of a clique and dissolves the board. We therefore find the core
        # as the largest clique (falling back to the longest path when no triangle
        # exists, as in a strictly nucleated cascade), and only then call a piece a
        # BRIDGE if it binds two or more core members from outside the core.
        self.contenders = self._contenders()
        benched = {c for c, _ in self.contenders}
        avail = {p for p in self.adj if p not in benched}

        self.backbone = self._core(avail, benched)
        core = set(self.backbone)

        self.bridges = {}
        self.pendants = {}
        for p in sorted(avail - core):
            hits = sorted(self.adj[p] & core, key=self.backbone.index)
            if len(hits) >= 2:
                self.bridges[p] = (hits[0], hits[1])
            elif len(hits) == 1:
                self.pendants[p] = hits[0]

    # ------------------------------------------------------------- topology
    def _core(self, avail, benched) -> List[str]:
        """The largest clique, or the longest path if there is no triangle.

        A SNARE bundle is a clique; a nucleated cascade is a chain. Both must work."""
        cliques: List[List[str]] = []

        def grow(clique, cands):
            if len(clique) >= 3:
                cliques.append(list(clique))
            for v in sorted(cands):
                grow(clique + [v], {u for u in cands if u > v and u in self.adj[v]})

        grow([], set(avail))
        if not cliques:
            return self._longest_path(avail, benched)

        # More than one clique may exist (the SNARE bundle and the complexin-bound
        # bundle are both triangles). The CORE is the largest, and among equals the
        # most tightly bound -- ranked by total interface contact coverage. A weakly
        # engaged accessory must not be mistaken for the core it decorates.
        def weight(c):
            return sum(self._cov(a, b)
                       for i, a in enumerate(c) for b in c[i + 1:])

        return max(cliques, key=lambda c: (len(c), weight(c)))

    def _longest_path(self, nodes, removed) -> List[str]:
        best: List[str] = []

        def walk(node, path):
            nonlocal best
            if len(path) > len(best):
                best = list(path)
            for n in sorted(self.adj[node] & nodes):
                if n not in path:
                    walk(n, path + [n])

        for start in sorted(nodes):
            walk(start, [start])
        return best

    def _bridges_free(self, benched) -> Dict[str, tuple]:
        """A BRIDGE binds two proteins that are themselves adjacent -- so it cannot
        be seated until the seam between them closes. Derived from the graph, not
        declared. Detected BEFORE the backbone, or it would be absorbed into it."""
        out = {}
        for p in sorted(self.adj):
            if p in benched:
                continue
            ns = sorted(self.adj[p] - benched)
            for a, b in itertools.combinations(ns, 2):
                if b in self.adj[a] and len(self.adj[a]) > 2 and len(self.adj[b]) > 2:
                    out[p] = (a, b)      # a and b are adjacent AND better connected
                    break
        return out

    def _contenders(self) -> List[tuple]:
        """Two partners sharing one site on X are NOT necessarily competing.

        If they are ADJACENT to one another, X lies in the groove BETWEEN them: X is
        a bridge, and both partners are present simultaneously. Complexin's central
        helix binds in the groove between two SNARE helices -- one interface, two
        partners, and that is exactly what makes it a bridge.

        They compete only if they are NOT adjacent, i.e. they cannot both be there.
        AP180's ANTH domain and SNAP25 both take VAMP2's SNARE motif, and AP180 does
        not bind SNAP25 -- so one excludes the other.
        """
        users: Dict[tuple, list] = defaultdict(list)
        for r in self.edges.itertuples():
            for me, other in ((r.protein_a, r.protein_b), (r.protein_b, r.protein_a)):
                st = self.site.get((me, other))
                if isinstance(st, str) and st:
                    users[(me, st)].append(other)

        out, seen = [], set()
        for (host, st), partners in users.items():
            if len(partners) < 2:
                continue
            # partners that are mutually adjacent coexist (a bridge sits between them)
            competing = [p for p in partners
                         if not any(q in self.adj[p] for q in partners if q != p)]
            if len(competing) < 2:
                continue
            winner = max(competing, key=lambda p: self.cov[frozenset((host, p))] or 0)
            for p in competing:
                if p != winner and p not in seen:
                    seen.add(p)
                    out.append((p, (host, st, winner)))
        return out

    def _style(self, n):
        f = self.meta.loc[n, "function"] if n in self.meta.index else "other"
        return PALETTE.get(f, PALETTE["other"])

    def _age(self, n):
        return AGE_RINGS.get(self.meta.loc[n, "age"], 1) if n in self.meta.index else 1

    def _cov(self, a, b):
        return self.cov.get(frozenset((a, b))) or 0

    # --------------------------------------------------------------- drawing
    def draw(self, out: str, dpi: int = 300):
        bb = self.backbone
        W, H, y0, x0 = 1.9, 1.55, 4.5, 1.2
        xs = {p: x0 + i * W for i, p in enumerate(bb)}
        yM = y0 + H / 2
        right = x0 + len(bb) * W

        fig, ax = plt.subplots(figsize=(max(11.0, 0.95 * right + 2.4), 7.8),
                               facecolor="white")
        ax.set_xlim(0.3, right + 1.2)
        ax.set_ylim(0.5, 8.6)
        ax.set_aspect("equal")
        ax.axis("off")

        SY = [y0 + H * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]

        def SX(p):
            return [xs[p] + W * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]

        for i, p in enumerate(bb):
            pc = Piece(xs[p], y0, W, H)
            if i > 0:
                pc.ladder("left", SY, set(range(self._cov(bb[i - 1], p))), SOCKET)
            if i < len(bb) - 1:
                pc.ladder("right", SY, set(range(self._cov(p, bb[i + 1]))), TAB)
            for b, (a, c) in self.bridges.items():
                if p == a:
                    pc.top(xs[p] + W * 0.72, TAB, KNOB)
                if p == c:
                    pc.top(xs[p] + W * 0.28, TAB, KNOB)
            kids = sorted(k for k, h in self.pendants.items() if h == p)
            for j, pd in enumerate(kids):
                lo = 0.10 + j * (0.80 / max(len(kids), 1))
                coords = [xs[p] + W * (lo + f * (0.80 / max(len(kids), 1)) * 0.8)
                          for f in (0.15, 0.35, 0.5, 0.65, 0.85)]
                pc.ladder("bottom", coords, set(range(self._cov(p, pd))), SOCKET)
            fc, ec = self._style(p)
            ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                 edgecolor=ec, lw=1.6, zorder=3))
            ax.text(xs[p] + W / 2, yM + 0.02, p, color=INK, fontsize=8.6,
                    ha="center", va="center", fontweight="bold", zorder=6)
            if i > 0:
                ax.text(xs[p] + W / 2, yM - 0.30,
                        f"{self._cov(bb[i - 1], p)}/{N_SUB}", color=ec,
                        fontsize=7.0, ha="center", fontweight="bold", zorder=6)
            _rings(ax, xs[p] + 0.26, y0 + H - 0.26, self._age(p))

        for b, (a, c) in self.bridges.items():
            if a not in xs or c not in xs:
                continue
            if bb.index(a) > bb.index(c):
                a, c = c, a
            bx = xs[a] + W * 0.50
            bw = (xs[c] - xs[a]) + W * 0.40
            pc = Piece(bx, y0 + H, bw, 0.92)
            pc.bottom(xs[a] + W * 0.72, SOCKET, KNOB)
            pc.bottom(xs[c] + W * 0.28, SOCKET, KNOB)
            fc, ec = self._style(b)
            ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                 edgecolor=ec, lw=1.5, zorder=2))
            ax.text(bx + bw / 2 + 0.18, y0 + H + 0.64, b, color=INK, fontsize=8.0,
                    ha="center", va="center", fontweight="bold", zorder=6)
            _rings(ax, bx + 0.24, y0 + H + 0.64, self._age(b))

        if self.bridges:
            ax.text(x0, y0 + H + 1.38, "BRIDGE PIECES \u2014 two sockets: unplaceable "
                    "until the seam beneath them closes", color="#3d7ab8",
                    fontsize=7.4, fontweight="bold")

        by_host: Dict[str, list] = defaultdict(list)
        for pd, host in self.pendants.items():
            by_host[host].append(pd)
        for host, pds in by_host.items():
            k = len(pds)
            pw = (W + 0.5) / k if k > 1 else W + 0.56
            for j, pd in enumerate(sorted(pds)):
                px = xs[host] - 0.25 + j * pw
                pc = Piece(px, y0 - 1.14, pw - 0.04, 1.08)
                # tabs must sit under the host's own sockets
                lo = 0.10 + j * (0.80 / k)
                coords = [xs[host] + W * (lo + f * (0.80 / k) * 0.8)
                          for f in (0.15, 0.35, 0.5, 0.65, 0.85)]
                pc.ladder("top", coords, set(range(self._cov(host, pd))), TAB)
                fc, ec = self._style(pd)
                ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                     edgecolor=ec, lw=1.5, zorder=2))
                ax.text(px + pw / 2 - 0.02, y0 - 0.52, pd, color=INK,
                        fontsize=7.0 if k > 1 else 7.8, ha="center", va="center",
                        fontweight="bold", zorder=6)
                ax.text(px + pw / 2 - 0.02, y0 - 0.82,
                        f"{self._cov(host, pd)}/{N_SUB}", color=ec, fontsize=6.8,
                        ha="center", fontweight="bold", zorder=6)
                _rings(ax, px + 0.22, y0 - 0.32, self._age(pd))

        if self.contenders:
            bw = max(len(self.contenders) * 2.6 + 0.5, 5.5)
            ax.add_patch(FancyBboxPatch((x0, 0.95), bw, 1.85,
                                        boxstyle="round,pad=0.07", facecolor="white",
                                        edgecolor=STOP, lw=1.3, ls=(0, (5, 3)), zorder=1))
            ax.text(x0 + 0.2, 2.52, "ALTERNATIVE OCCUPANCY \u2014 the same site, "
                    "another state", color=STOP, fontsize=7.6, fontweight="bold", zorder=6)
            ax.text(x0 + 0.2, 2.26, "overlap establishes NON-simultaneity; it does not "
                    "say whether that is competition or succession",
                    color=MUTED, fontsize=6.4, zorder=6)
            for i, (c, (host, site, occ)) in enumerate(self.contenders):
                cx = x0 + 0.25 + i * 2.6
                pc = Piece(cx, 1.15, 2.3, 0.85)
                pc.ladder("left", [1.15 + 0.85 * f
                                   for f in (0.2, 0.35, 0.5, 0.65, 0.8)],
                          set(range(self._cov(host, c))), TAB)
                fc, ec = self._style(c)
                ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                     edgecolor=ec, lw=1.3, zorder=3))
                ax.text(cx + 1.2, 1.74, c, color=INK, fontsize=8.0, ha="center",
                        fontweight="bold", zorder=6)
                ax.text(cx + 1.2, 1.42,
                        f"{self._cov(host, c)}/{N_SUB} \u00b7 \u2297 not with {occ}",
                        color=STOP, fontsize=6.4, ha="center", zorder=6)
                _rings(ax, cx + 0.24, 1.84, self._age(c))

        fig.tight_layout()
        fig.savefig(out, dpi=dpi, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out

    def report(self) -> str:
        L = [f"backbone : {' - '.join(self.backbone)}",
             "bridges  : " + (", ".join(f"{b} over {a}|{c}"
                                        for b, (a, c) in self.bridges.items()) or "none"),
             "pendants : " + (", ".join(f"{p} on {h}"
                                        for p, h in self.pendants.items()) or "none"),
             "benched  : " + (", ".join(f"{c} (site {s} on {h} held by {o})"
                                        for c, (h, s, o) in self.contenders) or "none"),
             "",
             "  A bridge binds two ADJACENT backbone proteins and cannot be seated",
             "  until the seam between them closes. That is the precedence constraint",
             "  counted in assembly.py -- derived here from the table, not declared."]
        return "\n".join(L)


def render(proteins_csv, interactions_csv, out, dpi=300):
    b = Board(pd.read_csv(proteins_csv), pd.read_csv(interactions_csv))
    b.draw(out, dpi=dpi)
    return b, out
