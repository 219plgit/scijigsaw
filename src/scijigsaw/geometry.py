"""Jigsaw geometry: complementary tabs and sockets.

A piece is a rectangle whose edges carry *features*: a tab (+1) protrudes, a
socket (-1) indents, and a flat edge carries nothing. Two pieces mate where a
tab meets a socket at the same physical point, which is why interfaces are
declared here in WORLD coordinates rather than as edge fractions: a tab and its
matching socket are then stated at the same location and cannot drift apart.

An interface may be drawn as a single feature or as a LADDER of N subsites, of
which a partner engages n. The ratio n/N is *interface contact coverage* -- a
visual proxy for the extent of interface engaged. It is not an affinity.
"""
from __future__ import annotations

import numpy as np

KNOB = 0.115      # radius of a whole-interface knob, as a fraction of edge length
SUB = 0.052       # radius of one subsite within a ladder

TAB, SOCKET, FLAT = +1, -1, 0


def edge_path(p0, p1, features=(), n_arc: int = 32) -> list:
    """Points along an edge from p0 to p1 carrying `features`.

    features: iterable of (s, kind[, knob]) with s in (0, 1) the position along
    the edge, kind in {TAB, SOCKET}. Traversal is counter-clockwise, so the
    outward normal is to the right of travel.
    """
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    v = p1 - p0
    L = float(np.hypot(*v))
    if L == 0:
        raise ValueError("degenerate edge")
    u = v / L
    normal = np.array([u[1], -u[0]])          # outward
    pts, cursor = [], 0.0
    for f in sorted(features, key=lambda g: g[0]):
        s, kind = f[0], f[1]
        knob = f[2] if len(f) > 2 else KNOB
        shape = f[3] if len(f) > 3 else "round"
        r = knob * L
        pts += [p0 + u * (L * t) for t in np.linspace(cursor, max(s - knob, cursor), 4)]
        c = p0 + u * (L * s)
        if shape == "round":
            # same sweep for tab and socket; `kind` flips the bulge direction, so
            # a tab (+1) is convex (out) and a socket (-1) is concave (a notch in)
            th = np.linspace(np.pi, 0, n_arc)
            pts += [c + u * (r * np.cos(t)) + normal * (kind * r * np.sin(t) * 1.30)
                    for t in th]
        elif shape == "wedge":
            # an isosceles tab/socket: a triangle keyed by size
            pts += [c - u * r,
                    c + normal * (kind * r * 1.45),
                    c + u * r]
        elif shape == "square":
            # a rectangular tab/socket with a slight neck, keyed by size
            d = normal * (kind * r * 1.30)
            pts += [c - u * r, c - u * (r * 0.62) + d,
                    c + u * (r * 0.62) + d, c + u * r]
        else:
            raise ValueError(f"unknown connector shape {shape!r}")
        cursor = s + knob
    pts += [p0 + u * (L * t) for t in np.linspace(cursor, 1.0, 4)]
    return pts


class Piece:
    """A jigsaw piece. Interfaces are declared in world coordinates."""

    __slots__ = ("x", "y", "w", "h", "edges")

    def __init__(self, x: float, y: float, w: float, h: float):
        self.x, self.y, self.w, self.h = x, y, w, h
        self.edges = {"bottom": [], "right": [], "top": [], "left": []}

    # -- one feature at a world coordinate --------------------------------
    def bottom(self, X, kind, knob=SUB, shape="round"):
        self.edges["bottom"].append(((X - self.x) / self.w, kind, knob, shape)); return self

    def right(self, Y, kind, knob=SUB, shape="round"):
        self.edges["right"].append(((Y - self.y) / self.h, kind, knob, shape)); return self

    def top(self, X, kind, knob=SUB, shape="round"):
        self.edges["top"].append(((self.x + self.w - X) / self.w, kind, knob, shape)); return self

    def left(self, Y, kind, knob=SUB, shape="round"):
        self.edges["left"].append(((self.y + self.h - Y) / self.h, kind, knob, shape)); return self

    # -- an interface as a ladder of subsites ------------------------------
    def ladder(self, edge: str, coords, engaged, kind, knob=SUB):
        """Place a feature at coords[i] for every i in `engaged`.

        The hub declares the FULL ladder as sockets; a partner declares tabs
        only at the subsites it engages, so unengaged subsites remain visibly
        empty. Coverage is then read by counting holes.
        """
        fn = {"bottom": self.bottom, "right": self.right,
              "top": self.top, "left": self.left}[edge]
        for i, c in enumerate(coords):
            if i in engaged:
                fn(c, kind, knob)
        return self

    def polygon(self) -> np.ndarray:
        x, y, w, h = self.x, self.y, self.w, self.h
        return np.array(
            edge_path((x, y), (x + w, y), self.edges["bottom"])
            + edge_path((x + w, y), (x + w, y + h), self.edges["right"])
            + edge_path((x + w, y + h), (x, y + h), self.edges["top"])
            + edge_path((x, y + h), (x, y), self.edges["left"]))


def coverage(engaged: int, n_sub: int) -> float:
    """Interface contact coverage, n/N.

    A visual proxy for the extent of interface engaged. NOT an affinity: a 5/5
    electrostatically mismatched interface may bind more weakly than a 2/5
    hydrophobic one. N is a display parameter (default 5).
    """
    if not 0 <= engaged <= n_sub:
        raise ValueError("engaged must lie in [0, n_sub]")
    return engaged / n_sub
