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
PANEL_LABELS = True   # set False to render the assembled board without A/B panel letters
SHOW_ALT_OCCUPANCY = True   # set False to render only the assembled complex (no alt-occupancy box)
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
        # n/N is a STRUCTURAL footprint fraction, supplied only for structure-derived
        # relations (extractor columns below). It is never derived from affinity; a
        # curated interactions.csv without these columns shows no n/N.
        have = "coverage" in e.columns
        e["coverage"] = e["coverage"] if have else None
        e["coverage_denom"] = (e["coverage_denom"] if "coverage_denom" in e.columns
                               else (N_SUB if have else None))
        self.edges = e

        self.adj: Dict[str, set] = defaultdict(set)
        self.cov: Dict[frozenset, int] = {}
        self.covN: Dict[frozenset, int] = {}
        self.site: Dict[tuple, str] = {}
        for r in e.itertuples():
            self.adj[r.protein_a].add(r.protein_b)
            self.adj[r.protein_b].add(r.protein_a)
            key = frozenset((r.protein_a, r.protein_b))
            cv = getattr(r, "coverage", None)
            if cv is not None and not (isinstance(cv, float) and math.isnan(cv)):
                self.cov[key] = int(cv)
                cn = getattr(r, "coverage_denom", None)
                self.covN[key] = int(cn) if cn and not (isinstance(cn, float) and math.isnan(cn)) else N_SUB
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
        # most functionally homogeneous -- the bundle whose members share a single
        # functional class (the fusion SNAREs), rather than a core decorated by a
        # regulator. This uses declared function metadata, not affinity.
        def weight(c):
            funcs = [self.meta.loc[v, "function"] for v in c if v in self.meta.index]
            return max((funcs.count(f) for f in set(funcs)), default=0)

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
            # the occupant is the most-integrated competitor (highest interface degree);
            # for the VAMP2 SNARE motif this is SNAP25. No coverage/affinity is used.
            winner = max(competing, key=lambda p: len(self.adj[p]))
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
        return self.cov.get(frozenset((a, b)))

    def _covN(self, a, b):
        return self.covN.get(frozenset((a, b)), N_SUB)

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
        if PANEL_LABELS:
            ax.text(0.5, 8.35, "A", fontsize=16, fontweight="bold", color=INK, va="top", zorder=10)

        SY = [y0 + H * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]

        def SX(p):
            return [xs[p] + W * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]

        _R = 0.22                      # absolute connector radius (clean single tabs/sockets)
        _kn = lambda edge: _R / edge   # -> knob as a fraction of that edge length

        # pendant layout: centred under each host, clear gaps, no cross-host overlap
        _PGAP = 0.16
        _byh: Dict[str, list] = defaultdict(list)
        for _pd, _host in self.pendants.items():
            _byh[_host].append(_pd)
        pend_layout: Dict[str, list] = {}
        for _host, _pds in _byh.items():
            _pds = sorted(_pds); _k = len(_pds)
            _span = W * 0.86
            _pw = (_span - (_k - 1) * _PGAP) / _k
            _sx = xs[_host] + (W - _span) / 2.0
            _rows = []
            for _j, _pd in enumerate(_pds):
                _px = _sx + _j * (_pw + _PGAP)
                _cx = _px + _pw / 2.0
                _half = _pw * 0.31
                _co = [_cx - _half + 2 * _half * _t for _t in (0.0, 0.25, 0.5, 0.75, 1.0)]
                _rows.append((_pd, _px, _pw, _co))
            pend_layout[_host] = _rows

        for i, p in enumerate(bb):
            pc = Piece(xs[p], y0, W, H)
            if i > 0:
                pc.left(yM, SOCKET, _kn(H))
            if i < len(bb) - 1:
                pc.right(yM, TAB, _kn(H))
            for b, (a, c) in self.bridges.items():
                if p == a:
                    pc.top(xs[p] + W * 0.72, TAB, _kn(W))
                if p == c:
                    pc.top(xs[p] + W * 0.28, TAB, _kn(W))
            for _pd, _px, _pw, _co in pend_layout.get(p, []):
                pc.bottom(_px + _pw / 2.0, TAB, _kn(W))
            fc, ec = self._style(p)
            ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                 edgecolor=ec, lw=1.6, zorder=3))
            ax.text(xs[p] + W / 2, yM + 0.02, p, color=INK, fontsize=8.6,
                    ha="center", va="center", fontweight="bold", zorder=6)
            if i > 0 and self._cov(bb[i - 1], p) is not None:
                ax.text(xs[p] + W / 2, yM - 0.30,
                        f"{self._cov(bb[i - 1], p)}/{self._covN(bb[i - 1], p)}", color=ec,
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
            pc.bottom(xs[a] + W * 0.72, SOCKET, _kn(bw))
            pc.bottom(xs[c] + W * 0.28, SOCKET, _kn(bw))
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

        for host, rows in pend_layout.items():
            single = len(rows) == 1
            for pd, px, pw, coords in rows:
                pc = Piece(px, y0 - 1.14, pw, 1.14)
                pc.top(px + pw / 2.0, SOCKET, _kn(pw))
                fc, ec = self._style(pd)
                ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                     edgecolor=ec, lw=1.5, zorder=2))
                cx = px + pw / 2.0
                ax.text(cx, y0 - 0.52, pd, color=INK,
                        fontsize=7.6 if single else 6.6, ha="center", va="center",
                        fontweight="bold", zorder=6)
                if self._cov(host, pd) is not None:
                    ax.text(cx, y0 - 0.82, f"{self._cov(host, pd)}/{self._covN(host, pd)}",
                            color=ec, fontsize=6.8, ha="center", fontweight="bold", zorder=6)
                _rings(ax, px + 0.20, y0 - 0.32, self._age(pd))

        if self.contenders and SHOW_ALT_OCCUPANCY:
            sites = {(h, s) for _, (h, s, o) in self.contenders}
            winners = {o for _, (h, s, o) in self.contenders}
            shared = len(sites) == 1 and len(winners) == 1
            if shared:
                h0, s0 = next(iter(sites)); w0 = next(iter(winners))
                sh = s0.replace("_", " ")
                chips = [(w0, True)] + [(c, False) for c, _ in self.contenders]
                bw = max(len(chips) * 2.05 + 0.5, 5.5)
                by, bh = 0.7, 2.35
                ax.add_patch(FancyBboxPatch((x0, by), bw, bh, boxstyle="round,pad=0.07",
                             facecolor="white", edgecolor=STOP, lw=1.3, ls=(0, (5, 3)), zorder=1))
                if PANEL_LABELS:
                    ax.text(0.5, by + bh, "B", fontsize=16, fontweight="bold", color=INK, va="top", zorder=10)
                ax.text(x0 + 0.2, by + bh - 0.24,
                        "ALTERNATIVE OCCUPANCY \u2014 one socket on %s\u2019s %s" % (h0, sh),
                        color=STOP, fontsize=7.3, fontweight="bold", zorder=6)
                ax.text(x0 + 0.2, by + bh - 0.50,
                        "%s fills it in the fusion complex; the others (marked \u2297) take the "
                        "same socket only in the retrieval state \u2014 one occupant at a time" % w0,
                        color=MUTED, fontsize=6.1, zorder=6)
                # the one shared socket, drawn as a labelled slot
                bx0, bx1 = x0 + 0.4, x0 + bw - 0.4
                sy = by + bh - 1.02
                ax.add_patch(FancyBboxPatch((bx0, sy - 0.17), bx1 - bx0, 0.34,
                             boxstyle="round,pad=0.02", facecolor=PALETTE["fusion"][0],
                             edgecolor=INK, lw=1.3, zorder=3))
                ax.text((bx0 + bx1) / 2, sy, "%s socket on %s" % (sh, h0), color=INK,
                        fontsize=6.7, ha="center", va="center", fontweight="bold", zorder=5)
                # candidates below, each with an arrow up into the one socket
                cyc = by + 0.5
                span = bx1 - bx0
                for i, (nm, win) in enumerate(chips):
                    cx = bx0 + span * (i + 0.5) / len(chips)
                    fc, ec = self._style(nm)
                    ax.add_patch(FancyBboxPatch((cx - 0.85, cyc - 0.29), 1.7, 0.58,
                                 boxstyle="round,pad=0.02", facecolor=fc, edgecolor=ec,
                                 lw=1.4 if win else 1.2,
                                 ls="solid" if win else (0, (3, 2)), zorder=3))
                    ax.text(cx, cyc + 0.04, nm, color=INK, fontsize=6.9, ha="center",
                            va="center", fontweight="bold", zorder=5)
                    ax.text(cx, cyc - 0.16, "fusion \u2013 occupies" if win else "retrieval",
                            color=(INK if win else STOP), fontsize=5.1, ha="center",
                            va="center", zorder=5)
                    ax.annotate("", xy=(cx, sy - 0.19), xytext=(cx, cyc + 0.31),
                                arrowprops=dict(arrowstyle="-|>", lw=1.2,
                                color=INK if win else STOP,
                                linestyle="solid" if win else "dashed"), zorder=4)
                    if not win:
                        ax.text(cx + 0.2, (sy - 0.19 + cyc + 0.31) / 2, "\u2297", color=STOP,
                                fontsize=8.5, ha="center", va="center", zorder=6)
            else:
                bw = max(len(self.contenders) * 2.6 + 0.5, 5.5)
                ax.add_patch(FancyBboxPatch((x0, 0.95), bw, 1.85, boxstyle="round,pad=0.07",
                             facecolor="white", edgecolor=STOP, lw=1.3, ls=(0, (5, 3)), zorder=1))
                ax.text(x0 + 0.2, 2.55, "ALTERNATIVE OCCUPANCY \u2014 a site already held, "
                        "another state", color=STOP, fontsize=7.3, fontweight="bold", zorder=6)
                ax.text(x0 + 0.2, 2.30, "each piece binds a site an occupant already holds \u2014 "
                        "only one at a time", color=MUTED, fontsize=6.3, zorder=6)
                for i, (c, (host, site, occ)) in enumerate(self.contenders):
                    cx = x0 + 0.25 + i * 2.6
                    pc = Piece(cx, 1.10, 2.3, 0.82)
                    pc.left(1.10 + 0.82 / 2.0, TAB, _kn(0.82))
                    fc, ec = self._style(c)
                    ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc,
                                         edgecolor=ec, lw=1.3, zorder=3))
                    ax.text(cx + 1.28, 1.72, c, color=INK, fontsize=8.0, ha="center",
                            fontweight="bold", zorder=6)
                    ax.text(cx + 1.28, 1.47, "binds %s \u00b7 %s" % (host, site.replace("_", " ")),
                            color=MUTED, fontsize=5.9, ha="center", zorder=6)
                    ax.text(cx + 1.28, 1.25, "\u2297 excludes %s" % occ, color=STOP,
                            fontsize=6.2, ha="center", zorder=6)
                    _rings(ax, cx + 0.24, 1.82, self._age(c))


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
