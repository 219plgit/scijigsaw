"""Board renderer: an interaction table in, a jigsaw board out.

This is the second tier of the architecture. `extract.py` derives interfaces from
structures and writes an interaction table; this module lays that table out as a
board and draws it. The two are decoupled deliberately: a curated table renders
without any structure at all.

    scijigsaw-render examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
        --out vamp2_board.svg

WHAT THE GEOMETRY ENCODES
    shape     tabs and sockets: a piece mates only where an interface is mapped
    colour    functional class
    rings     evolutionary age (more rings = older)
    n/N       interface contact coverage -- subsites engaged, NOT affinity

WHAT IT REFUSES
    A partner with no site assignment is drawn dashed and flagged: a private site
    is ASSUMED, not measured. Partners sharing a site cannot be seated together;
    the board shows one admissible occupancy and the alternatives sit beside it,
    because overlap establishes non-simultaneity, not competition (an enzyme may
    drive succession instead).
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

from .geometry import SOCKET, TAB, KNOB, SUB, Piece

# ----------------------------------------------------------------- appearance
INK, MUTED, RULE = "#1a1d24", "#5b6675", "#c7cdd6"
STOP = "#c0392b"
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
# Kd (nM) -> subsites engaged.  A DECLARED CONVENTION, not a law of nature:
# contact coverage proxies buried interface, not affinity.
KD_LADDER = [(1.0, 5), (100.0, 4), (1_000.0, 3), (100_000.0, 2), (math.inf, 1)]


def kd_to_coverage(kd_nM) -> int | None:
    if kd_nM is None or (isinstance(kd_nM, float) and math.isnan(kd_nM)):
        return None
    for limit, n in KD_LADDER:
        if float(kd_nM) < limit:
            return n
    return 1


def _rings(ax, x, y, n, colour=INK, r0=0.042, z=9):
    ax.add_patch(plt.Circle((x, y), r0 * (n + 1), facecolor="white",
                            edgecolor=colour, lw=0.7, zorder=z))
    for k in range(1, n + 1):
        ax.add_patch(plt.Circle((x, y), r0 * k, facecolor="none",
                                edgecolor=colour, lw=0.8, zorder=z + 1))


class Board:
    """Lay out and draw a hub-and-partner board from a table."""

    def __init__(self, proteins: pd.DataFrame, interactions: pd.DataFrame):
        self.meta = proteins.set_index("name")
        self.edges = interactions.copy()
        if "kd_nM" in self.edges:
            self.edges["coverage"] = self.edges["kd_nM"].apply(kd_to_coverage)
        else:
            self.edges["coverage"] = None
        for col in ("site_on_a", "site_on_b"):
            if col not in self.edges:
                self.edges[col] = np.nan
        self.edges["assumed_site"] = (self.edges["site_on_a"].isna()
                                      | (self.edges["site_on_a"] == ""))
        self.hub = self._pick_hub()
        self.sites = self._sites_on_hub()

    # ---------------------------------------------------------------- model
    def _pick_hub(self) -> str:
        deg: Dict[str, int] = defaultdict(int)
        for r in self.edges.itertuples():
            deg[r.protein_a] += 1
            deg[r.protein_b] += 1
        return max(deg, key=deg.get)

    def _sites_on_hub(self) -> Dict[str, list]:
        sites: Dict[str, list] = defaultdict(list)
        for r in self.edges.itertuples():
            for me, other, site in ((r.protein_a, r.protein_b, r.site_on_a),
                                    (r.protein_b, r.protein_a, r.site_on_b)):
                if me != self.hub:
                    continue
                key = site if isinstance(site, str) and site else f"__assumed__{other}"
                sites[key].append((other, r.coverage, bool(r.assumed_site)))
        return dict(sites)

    def contested(self) -> Dict[str, list]:
        return {s: v for s, v in self.sites.items() if len(v) > 1}

    def _style(self, name):
        fn = self.meta.loc[name, "function"] if name in self.meta.index else "other"
        return PALETTE.get(fn, PALETTE["other"])

    def _age(self, name):
        if name not in self.meta.index:
            return 1
        return AGE_RINGS.get(self.meta.loc[name, "age"], 1)

    # ---------------------------------------------------------------- render
    def draw(self, out: str, title: str | None = None, dpi: int = 300):
        fig, ax = plt.subplots(figsize=(11, 7), facecolor="white")
        ax.set_xlim(0, 16); ax.set_ylim(0, 10)
        ax.set_aspect("equal"); ax.axis("off")

        HX, HY, HW, HH = 6.4, 4.6, 3.0, 2.0
        yM = HY + HH / 2

        # one hub edge per site (up to four)
        slots = ["left", "right", "top", "bottom"]
        site_slot = {s: slots[i % 4] for i, s in enumerate(self.sites)}

        def ladder_coords(edge):
            if edge in ("left", "right"):
                return [HY + HH * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]
            return [HX + HW * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]

        hub = Piece(HX, HY, HW, HH)
        for s, slot in site_slot.items():
            hub.ladder(slot, ladder_coords(slot), set(range(N_SUB)), SOCKET)
        fc, ec = self._style(self.hub)
        ax.add_patch(Polygon(hub.polygon(), closed=True, facecolor=fc,
                             edgecolor=ec, lw=1.8, zorder=3))
        ax.text(HX + HW / 2, yM, self.hub, color=INK, fontsize=13, ha="center",
                va="center", fontweight="bold", zorder=6)
        _rings(ax, HX + 0.30, HY + HH - 0.30, self._age(self.hub))

        PW, PH = 2.8, 1.9                    # partner size
        # partners must be FLUSH with the hub: a tab has to enter a socket.
        anchor = {"left":   (HX - PW,      HY + (HH - PH) / 2, "right"),
                  "right":  (HX + HW,      HY + (HH - PH) / 2, "left"),
                  "top":    (HX + (HW - PW) / 2, HY + HH,      "bottom"),
                  "bottom": (HX + (HW - PW) / 2, HY - PH,      "top")}

        bench = []
        for s, users in self.sites.items():
            seated, *others = sorted(users, key=lambda t: -(t[1] or 0))
            name, cov, assumed = seated
            px, py, pedge = anchor[site_slot[s]]
            p = Piece(px, py, PW, PH)
            coords = (ladder_coords("left") if pedge in ("left", "right")
                      else ladder_coords("top"))
            # partner tabs must sit at the hub's world coordinates
            # the partner's tabs must sit at the HUB's socket coordinates
            coords = ladder_coords(site_slot[s])
            p.ladder(pedge, coords, set(range(cov or 0)), TAB)
            f2, e2 = self._style(name)
            ax.add_patch(Polygon(p.polygon(), closed=True, facecolor=f2,
                                 edgecolor=e2, lw=1.5,
                                 linestyle=(0, (5, 3)) if assumed else "-",
                                 zorder=2))
            ax.text(px + PW / 2, py + PH * 0.72, name, color=INK, fontsize=10,
                    ha="center", fontweight="bold", zorder=6)
            ax.text(px + PW / 2, py + PH * 0.47, f"{cov or 0}/{N_SUB}", color=e2,
                    fontsize=8.4, ha="center", fontweight="bold", zorder=6)
            label = "site assumed" if assumed else s
            ax.text(px + PW / 2, py + PH * 0.28, label, color=MUTED, fontsize=6.6,
                    ha="center", style="italic", zorder=6)
            _rings(ax, px + 0.30, py + PH - 0.28, self._age(name))
            for on, oc, oa in others:
                bench.append((on, oc, s, name))

        # alternatives: same site, other occupancy (NOT losers)
        if bench:
            ax.add_patch(FancyBboxPatch((0.5, 0.35), 15.0, 1.6,
                                        boxstyle="round,pad=0.08",
                                        facecolor="white", edgecolor=STOP,
                                        lw=1.3, linestyle=(0, (5, 3)), zorder=1))
            ax.text(0.8, 1.66, "ALTERNATIVE OCCUPANCY \u2014 the same site, another state",
                    color=STOP, fontsize=8.6, fontweight="bold", zorder=6)
            ax.text(0.8, 1.40, "overlap establishes that these cannot be seated "
                    "SIMULTANEOUSLY; it does not say whether that is competition "
                    "or temporal succession.", color=MUTED, fontsize=6.8, zorder=6)
            for i, (name, cov, s, occ) in enumerate(bench[:4]):
                bx = 0.9 + i * 3.7
                q = Piece(bx, 0.55, 2.6, 0.72)
                q.ladder("left", [0.55 + 0.72 * f for f in
                                  (0.20, 0.35, 0.50, 0.65, 0.80)],
                         set(range(cov or 0)), TAB)
                f2, e2 = self._style(name)
                ax.add_patch(Polygon(q.polygon(), closed=True, facecolor=f2,
                                     edgecolor=e2, lw=1.2, alpha=0.7, zorder=3))
                ax.text(bx + 1.35, 1.02, name, color=INK, fontsize=8.4,
                        ha="center", fontweight="bold", zorder=6)
                ax.text(bx + 1.35, 0.72, f"{cov or 0}/{N_SUB} \u00b7 \u2297 not with {occ}",
                        color=STOP, fontsize=6.4, ha="center", zorder=6)

        if title:
            ax.text(0.5, 9.5, title, color=INK, fontsize=14, fontweight="bold")

        fig.tight_layout()
        fig.savefig(out, dpi=dpi, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out

    # ------------------------------------------------------------- reporting
    def report(self) -> str:
        lines = [f"hub: {self.hub}", ""]
        lines.append("sites on the hub")
        for s, users in self.sites.items():
            tag = "   <-- CONTESTED" if len(users) > 1 else ""
            shown = s if not s.startswith("__assumed__") else "(site assumed)"
            lines.append("  " + f"{shown:26}" +
                         ", ".join(f"{u} {c or 0}/{N_SUB}" for u, c, _ in users) + tag)
        con = self.contested()
        if con:
            lines += ["", "not simultaneous (same site):"]
            for s, users in con.items():
                lines.append("  " + "  XOR  ".join(u for u, _, _ in users))
            lines += ["", "  Overlap says these cannot bind AT THE SAME TIME.",
                      "  It does NOT say why: competition and temporal succession",
                      "  give identical geometry."]
        assumed = self.edges[self.edges.assumed_site]
        if len(assumed):
            lines += ["", f"[!] {len(assumed)} edge(s) have no site annotation.",
                      "    A private site is ASSUMED and drawn dashed."]
        return "\n".join(lines)


def render(proteins_csv: str, interactions_csv: str, out: str,
           title: str | None = None, dpi: int = 300):
    """Render a board from an interaction table. Returns (Board, path)."""
    b = Board(pd.read_csv(proteins_csv), pd.read_csv(interactions_csv))
    b.draw(out, title=title, dpi=dpi)
    return b, out
