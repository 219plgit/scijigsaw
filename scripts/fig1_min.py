#!/usr/bin/env python3
"""Figure 1 (minimal) -- rule, channels, pipeline. Labels only; the argument
lives in the caption."""
import os, sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType (no Type 3)
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'
os.makedirs(OUT, exist_ok=True)
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

INK, MUTED, RULE = "#1a1d24", "#5b6675", "#c7cdd6"
CORE, CORE_ED = "#d7dbe2", "#78828f"
SPEC, SPEC_ED = "#cfe0f2", "#3d7ab8"
RETR, RETR_ED = "#cbe9e3", "#1f7a6f"
SYN, SYN_ED = "#fde9b0", "#b8860b"
STOP_ED = "#c0392b"
PANEL = "#f7f8fa"
KNOB = 0.115


def _edge(p0, p1, feats=(), n=32):
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


def pc(x, y, w, h, b=(), r=(), t=(), l=()):
    return np.array(_edge((x, y), (x + w, y), b) + _edge((x + w, y), (x + w, y + h), r)
                    + _edge((x + w, y + h), (x, y + h), t) + _edge((x, y + h), (x, y), l))


def rings(ax, x, y, n, col=INK, r0=0.048, z=9):
    ax.add_patch(plt.Circle((x, y), r0 * (n + 1), facecolor="white", edgecolor=col,
                            lw=0.7, zorder=z))
    for k in range(1, n + 1):
        ax.add_patch(plt.Circle((x, y), r0 * k, facecolor="none", edgecolor=col,
                                lw=0.8, zorder=z + 1))


fig, ax = plt.subplots(figsize=(10, 5.2), facecolor="white")
ax.set_xlim(0, 16); ax.set_ylim(2.6, 11); ax.set_aspect("equal"); ax.axis("off")

# ---------------- A : the rule ----------------------------------------------
ax.add_patch(FancyBboxPatch((0.35, 7.10), 6.4, 3.35, boxstyle="round,pad=0.05",
                            facecolor=PANEL, edgecolor=RULE, lw=1.0, zorder=1))
ax.text(0.62, 10.05, "A", color=INK, fontsize=11, fontweight="bold", zorder=6)

ax.add_patch(Polygon(pc(0.95, 8.15, 1.55, 1.30, r=[(0.5, -1)]), closed=True,
                     facecolor=CORE, edgecolor=CORE_ED, lw=1.4, zorder=3))
ax.add_patch(Polygon(pc(2.55, 8.15, 1.45, 1.30, l=[(0.5, +1)]), closed=True,
                     facecolor=SPEC, edgecolor=SPEC_ED, lw=1.4, zorder=3))
ax.text(2.48, 7.70, "\u2713", color=SPEC_ED, fontsize=15, ha="center", zorder=6)

ax.add_patch(Polygon(pc(4.75, 8.15, 1.45, 1.30), closed=True, facecolor="#e6e6e6",
                     edgecolor="#9aa0a6", lw=1.4, zorder=3))
ax.text(5.48, 7.70, "\u2717", color="#9aa0a6", fontsize=15, ha="center", zorder=6)

# ---------------- B : four channels -----------------------------------------
ax.add_patch(FancyBboxPatch((7.15, 7.10), 8.5, 3.35, boxstyle="round,pad=0.05",
                            facecolor=PANEL, edgecolor=RULE, lw=1.0, zorder=1))
ax.text(7.42, 10.05, "B", color=INK, fontsize=11, fontweight="bold", zorder=6)

ax.add_patch(Polygon(pc(8.15, 8.30, 1.15, 1.00, r=[(0.5, -1)]), closed=True,
                     facecolor=CORE, edgecolor=CORE_ED, lw=1.3, zorder=3))
ax.text(8.72, 7.85, "shape", color=INK, fontsize=8.4, ha="center",
        fontweight="bold", zorder=6)

for i, (fc, ec) in enumerate([(CORE, CORE_ED), (SPEC, SPEC_ED), (RETR, RETR_ED),
                              (SYN, SYN_ED)]):
    ax.add_patch(Rectangle((10.05 + i * 0.42, 8.55), 0.34, 0.50, facecolor=fc,
                           edgecolor=ec, lw=1.1, zorder=3))
ax.text(10.90, 7.85, "colour", color=INK, fontsize=8.4, ha="center",
        fontweight="bold", zorder=6)

for i, n in enumerate([3, 2, 1]):
    rings(ax, 12.55 + i * 0.52, 8.80, n)
ax.text(13.07, 7.85, "rings", color=INK, fontsize=8.4, ha="center",
        fontweight="bold", zorder=6)

ax.add_patch(Polygon(pc(14.45, 8.30, 0.95, 1.00,
                        r=[(0.28, -1, 0.14), (0.52, -1, 0.14), (0.76, -1, 0.14)]),
                     closed=True, facecolor=SYN, edgecolor=SYN_ED, lw=1.3, zorder=3))
ax.text(14.92, 7.85, "n / N", color=INK, fontsize=8.4, ha="center",
        fontweight="bold", zorder=6)

# ---------------- C : pipeline ----------------------------------------------
ax.add_patch(FancyBboxPatch((0.35, 3.00), 15.3, 3.55, boxstyle="round,pad=0.05",
                            facecolor=PANEL, edgecolor=RULE, lw=1.0, zorder=1))
ax.text(0.62, 6.15, "C", color=INK, fontsize=11, fontweight="bold", zorder=6)

ax.add_patch(FancyBboxPatch((0.85, 3.55), 6.6, 2.05, boxstyle="round,pad=0.05",
                            facecolor="white", edgecolor=RETR_ED, lw=1.2,
                            ls=(0, (5, 3)), zorder=2))
ax.text(1.15, 5.30, "structure-derived", color=RETR_ED, fontsize=6.6, fontweight="bold",
        zorder=6)
for x, w, lab, fc in [(1.15, 1.60, "structures", "white"),
                      (3.25, 2.35, "jigsaw_extract", RETR),
                      (6.05, 1.10, "sites", "white")]:
    ax.add_patch(FancyBboxPatch((x, 3.80), w, 0.95, boxstyle="round,pad=0.04",
                                facecolor=fc, edgecolor=RETR_ED, lw=1.1, zorder=3))
    ax.text(x + w / 2, 4.27, lab, color=INK, fontsize=7.8, ha="center",
            va="center", fontweight="bold", zorder=6)
for a, b in [(2.80, 3.20), (5.65, 6.00)]:
    ax.add_patch(FancyArrowPatch((a, 4.27), (b, 4.27), arrowstyle="-|>",
                                 mutation_scale=11, lw=1.2, color=RETR_ED, zorder=5))

ax.add_patch(FancyArrowPatch((7.55, 4.27), (8.75, 4.27), arrowstyle="-|>",
                             mutation_scale=14, lw=1.7, color=INK, zorder=5))

ax.add_patch(FancyBboxPatch((8.90, 3.55), 6.4, 2.05, boxstyle="round,pad=0.05",
                            facecolor="white", edgecolor=SPEC_ED, lw=1.2,
                            ls=(0, (5, 3)), zorder=2))
ax.text(9.20, 5.30, "curated-table", color=SPEC_ED, fontsize=6.6, fontweight="bold",
        zorder=6)
for x, w, lab, fc in [(9.20, 1.85, "interactions", "white"),
                      (11.45, 2.05, "jigsaw", SPEC),
                      (13.90, 1.25, "puzzle", "white")]:
    ax.add_patch(FancyBboxPatch((x, 3.80), w, 0.95, boxstyle="round,pad=0.04",
                                facecolor=fc, edgecolor=SPEC_ED, lw=1.1, zorder=3))
    ax.text(x + w / 2, 4.27, lab, color=INK, fontsize=7.8, ha="center",
            va="center", fontweight="bold", zorder=6)
for a, b in [(11.05, 11.40), (13.50, 13.85)]:
    ax.add_patch(FancyArrowPatch((a, 4.27), (b, 4.27), arrowstyle="-|>",
                                 mutation_scale=11, lw=1.2, color=SPEC_ED, zorder=5))

fig.savefig(os.path.join(OUT, "Figure1_tool.png"), dpi=300, facecolor="white", bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure1_tool.pdf"), facecolor="white", bbox_inches="tight")
print(f"  wrote Figure1_tool.pdf / .png")
