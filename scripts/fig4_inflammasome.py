#!/usr/bin/env python3
"""Figure 4 -- the NLRP3 inflammasome board.

The point of the figure is the SHAPE. The VAMP2 board (Fig. 2) is a HUB with
spokes: shallow, poset depth 3, pruned 15x. The inflammasome is a CHAIN: each
piece can only be placed once the one above it is present, poset depth 9, and
just TWO of 3,628,800 orders survive.

Labels only; the argument lives in the caption.

Assembly (Lu et al., Cell 2014; Xiao et al., Annu Rev Immunol 2023):
  NLRP3 -> NEK7 -> inflammasome disk -> NLRP3-PYD filament -> nucleates
  ASC-PYD filament -> ASC-CARD filament -> nucleates caspase-1 CARD filament
  -> caspase domain dimerises + autoprocesses -> cleaves GSDMD and pro-IL-1beta
"""
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'
os.makedirs(OUT, exist_ok=True)
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

INK, MUTED = "#1a1d24", "#5b6675"
SENSOR, SENSOR_ED = "#d7dbe2", "#78828f"      # NLRP3 / NEK7 / disk
ADAPT, ADAPT_ED = "#cfe0f2", "#3d7ab8"        # ASC
EFFECT, EFFECT_ED = "#cbe9e3", "#1f7a6f"      # caspase-1
SUBST, SUBST_ED = "#fbdc8a", "#a8760a"        # GSDMD / IL-1beta
KNOB = 0.115


def _edge(p0, p1, feats=(), n=30):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    v = p1 - p0
    L = np.hypot(*v)
    u = v / L
    nrm = np.array([u[1], -u[0]])
    pts, cur = [], 0.0
    for f in sorted(feats):
        s, tab = f[0], f[1]
        kn = f[2] if len(f) > 2 else KNOB
        r = kn * L
        pts += [p0 + u * (L * t) for t in np.linspace(cur, max(s - kn, cur), 4)]
        c = p0 + u * (L * s)
        th = (np.linspace(np.pi, 0, n) if tab > 0
              else np.linspace(np.pi, 2 * np.pi, n))
        pts += [c + u * (r * np.cos(t)) + nrm * (tab * r * np.sin(t) * 1.3) for t in th]
        cur = s + kn
    pts += [p0 + u * (L * t) for t in np.linspace(cur, 1.0, 4)]
    return pts


def piece(x, y, w, h, top=None, bottom=None):
    """top/bottom: +1 tab, -1 socket, None flat. Placed at the mid-point."""
    b = [(0.5, bottom)] if bottom else []
    t = [(0.5, top)] if top else []
    return np.array(_edge((x, y), (x + w, y), b)
                    + _edge((x + w, y), (x + w, y + h), [])
                    + _edge((x + w, y + h), (x, y + h), t)
                    + _edge((x, y + h), (x, y), []))


fig, ax = plt.subplots(figsize=(7.4, 9.4), facecolor="white")
ax.set_xlim(0.2, 7.6); ax.set_ylim(0.2, 10.6)
ax.set_aspect("equal"); ax.axis("off")

W, H, X = 3.05, 0.78, 1.35
GAP = 0.10

chain = [
    ("NLRP3",                  "sensor",                       SENSOR, SENSOR_ED),
    ("NEK7",                   "licensing kinase",             SENSOR, SENSOR_ED),
    ("NLRP3 disk",             "oligomer",                     SENSOR, SENSOR_ED),
    ("NLRP3 PYD filament",     "nucleates \u2193",             SENSOR, SENSOR_ED),
    ("ASC PYD filament",       "adaptor",                      ADAPT,  ADAPT_ED),
    ("ASC CARD filament",      "nucleates \u2193",             ADAPT,  ADAPT_ED),
    ("Caspase-1 CARD filament", "effector",                    EFFECT, EFFECT_ED),
    ("Caspase-1 catalytic",    "dimerises \u00b7 autocleaves", EFFECT, EFFECT_ED),
]

y = 10.05
for i, (name, sub, fc, ec) in enumerate(chain):
    y -= H + GAP
    top = -1 if i > 0 else None          # socket receiving the piece above
    bot = +1                             # tab into the piece below
    ax.add_patch(Polygon(piece(X, y, W, H, top=top, bottom=bot), closed=True,
                         facecolor=fc, edgecolor=ec, lw=1.5, zorder=3))
    ax.text(X + W / 2, y + H / 2 + 0.05, name, color=INK, fontsize=8.8,
            ha="center", va="center", fontweight="bold", zorder=6)

# ---- the fork: two substrates, both requiring the catalytic domain ----------
yb = y - H - GAP - 0.28
for k, (nm, sb) in enumerate([("GSDMD", "pore \u00b7 pyroptosis"),
                              ("pro-IL-1\u03b2", "cytokine")]):
    xb = X - 0.55 + k * 2.05
    ax.add_patch(Polygon(piece(xb, yb, 1.60, H, top=-1), closed=True,
                         facecolor=SUBST, edgecolor=SUBST_ED, lw=1.4, zorder=3))
    ax.text(xb + 0.80, yb + H / 2 + 0.04, nm, color="#3a2c00", fontsize=8.4,
            ha="center", va="center", fontweight="bold", zorder=6)
for k in range(2):
    xb = X - 0.55 + k * 2.05 + 0.80
    ax.add_patch(FancyArrowPatch((X + W / 2, y - 0.02), (xb, yb + H + 0.16),
                                 connectionstyle="arc3,rad=%0.2f" % (0.18 if k else -0.18),
                                 arrowstyle="-", lw=1.1, color=SUBST_ED,
                                 ls=(0, (3, 2)), zorder=2))

# ---- the contrast, stated as geometry ---------------------------------------
ax.annotate("", xy=(0.95, 1.30), xytext=(0.95, 9.95),
            arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#9aa4b4",
                            mutation_scale=14))
ax.text(0.62, 5.6, "each piece can only be placed\nonce the one above it is present",
        color=MUTED, fontsize=7.4, rotation=90, ha="center", va="center",
        linespacing=1.5)

ax.add_patch(FancyBboxPatch((4.75, 8.15), 2.65, 1.85, boxstyle="round,pad=0.07",
                            facecolor="#f6f7f9", edgecolor="#c7cdd6", lw=1.0,
                            zorder=1))
ax.text(6.08, 9.72, "a CHAIN, not a hub", color=INK, fontsize=8.8,
        ha="center", fontweight="bold", zorder=6)
for j, (k, v) in enumerate([("subunits", "10"), ("poset depth", "9"),
                            ("orders n!", "3,628,800"),
                            ("permitted", "2")]):
    yy = 9.38 - j * 0.30
    ax.text(4.92, yy, k, color=MUTED, fontsize=7.4, va="center", zorder=6)
    ax.text(7.25, yy, v, color=INK, fontsize=7.6, va="center", ha="right",
            fontweight="bold", zorder=6)

fig.savefig(os.path.join(OUT, "FigureS1_inflammasome.png"), dpi=300, facecolor="white", bbox_inches="tight")
fig.savefig(os.path.join(OUT, "FigureS1_inflammasome.pdf"), facecolor="white", bbox_inches="tight")
print(f"  wrote FigureS1_inflammasome.pdf / .png")
