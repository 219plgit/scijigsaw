#!/usr/bin/env python3
"""Figure 2 -- the VAMP2 board. Labels only; the argument lives in the caption.

Completeness (asked and answered before submission):
  IN   VAMP2, SNAP25, Syntaxin-1A     core fusion, ancient
       Munc18-1                       binds Syntaxin ONLY. Ancient (Sec1p) yet
                                      functionally 'specialisation' -> the exemplar
                                      of colour and rings disagreeing
       Complexin, Synaptotagmin-1     SPLIT into two bridge pieces (distinct
                                      proteins, distinct interfaces)
       Synaptophysin                  VAMP2's THIRD socket (vesicle face)
       AP180, CALM                    contest the SNARE motif -> benched
       SNCA                           VAMP2's N-terminal socket, 2/5
  OUT  SNCB/SNCG (paralogues -> caption); MOSPD/VAP (separate board, Fig 1);
       Rab3A/RIM/Munc13 (not direct VAMP2 partners)
  CANNOT BE DRAWN
       NSF / alpha-SNAP disassembles the cis-SNARE complex and thereby VACATES the
       contested socket -- it is what makes 'fusion XOR retrieval' switchable. A
       jigsaw encodes binding, competition and assembly order; it cannot encode
       CATALYSIS, i.e. an enzyme whose function is to remove another piece. This is
       stated as a limit of the representation rather than patched with ad-hoc
       notation.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

INK, MUTED = "#1a1d24", "#5b6675"
CORE, CORE_ED = "#d7dbe2", "#78828f"
SPEC, SPEC_ED = "#cfe0f2", "#3d7ab8"
RETR, RETR_ED = "#cbe9e3", "#1f7a6f"
SYN, SYN_ED = "#fbdc8a", "#a8760a"
STOP_ED = "#c0392b"
KNOB, SUB = 0.115, 0.052


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


class P:
    def __init__(s, x, y, w, h):
        s.x, s.y, s.w, s.h = x, y, w, h
        s.e = {"b": [], "r": [], "t": [], "l": []}

    def add(s, k, f, tab, kn=SUB):
        s.e[k].append((f, tab, kn)); return s

    def bottom(s, X, t, kn=SUB): return s.add("b", (X - s.x) / s.w, t, kn)
    def right(s, Y, t, kn=SUB): return s.add("r", (Y - s.y) / s.h, t, kn)
    def top(s, X, t, kn=SUB): return s.add("t", (s.x + s.w - X) / s.w, t, kn)
    def left(s, Y, t, kn=SUB): return s.add("l", (s.y + s.h - Y) / s.h, t, kn)

    def ladder(s, e, coords, eng, tab):
        f = {"b": s.bottom, "r": s.right, "t": s.top, "l": s.left}[e]
        for i, c in enumerate(coords):
            if i in eng:
                f(c, tab)
        return s

    def poly(s):
        x, y, w, h = s.x, s.y, s.w, s.h
        return np.array(_edge((x, y), (x + w, y), s.e["b"])
                        + _edge((x + w, y), (x + w, y + h), s.e["r"])
                        + _edge((x + w, y + h), (x, y + h), s.e["t"])
                        + _edge((x, y + h), (x, y), s.e["l"]))

    def draw(s, ax, fc, ec, z=3, lw=1.5):
        ax.add_patch(Polygon(s.poly(), closed=True, facecolor=fc, edgecolor=ec,
                             lw=lw, zorder=z)); return s


def rings(ax, x, y, n, col=INK, r0=0.042, z=9):
    ax.add_patch(plt.Circle((x, y), r0 * (n + 1), facecolor="white", edgecolor=col,
                            lw=0.7, zorder=z))
    for k in range(1, n + 1):
        ax.add_patch(plt.Circle((x, y), r0 * k, facecolor="none", edgecolor=col,
                                lw=0.8, zorder=z + 1))


fig, ax = plt.subplots(figsize=(10.5, 7.0), facecolor="white")
ax.set_xlim(0.75, 10.5); ax.set_ylim(1.0, 8.1); ax.set_aspect("equal"); ax.axis("off")

Y, H, W = 4.35, 1.55, 1.80
xV, xS, xQ = 3.05, 4.85, 6.65
yM = Y + H / 2

SN = [Y + H * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]
NT = [Y + H * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)]
SP = [xV + W * f for f in (0.25, 0.40, 0.55, 0.70, 0.85)]
MU = [Y + H * f for f in (0.24, 0.44, 0.64, 0.84)]
ALL5, ANTH3, SYN2, SYP3 = {0, 1, 2, 3, 4}, {0, 1, 2}, {2, 3}, {1, 2, 3}
cpxV, cpxS, sytS, sytQ = 3.95, 5.30, 5.90, 7.25

# ---------------- CORE -------------------------------------------------------
v = P(xV, Y, W, H).top(cpxV, +1, KNOB)
v.ladder("r", SN, ALL5, -1)
v.ladder("l", NT, ALL5, -1)
v.ladder("b", SP, ALL5, -1)
v.draw(ax, CORE, CORE_ED)
ax.text(xV + W / 2 - 0.12, yM, "VAMP2", color=INK, fontsize=9.0, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, xV + 0.26, Y + H - 0.25, 3)

s = P(xS, Y, W, H).right(yM, +1, KNOB).top(cpxS, +1, KNOB).top(sytS, +1, KNOB)
s.ladder("l", SN, ALL5, +1)
s.draw(ax, CORE, CORE_ED)
ax.text(xS + W / 2 + 0.12, yM + 0.16, "SNAP25", color=INK, fontsize=8.6,
        ha="center", va="center", fontweight="bold", zorder=6)
ax.text(xS + W / 2 + 0.12, yM - 0.26, "5/5", color=CORE_ED, fontsize=7.6,
        ha="center", va="center", fontweight="bold", zorder=6)
rings(ax, xS + 0.26, Y + H - 0.25, 3)

q = P(xQ, Y, W, H).left(yM, -1, KNOB).top(sytQ, +1, KNOB)
q.ladder("r", MU, {0, 1, 2, 3}, +1)
q.draw(ax, CORE, CORE_ED)
ax.text(xQ + W / 2 + 0.08, yM, "Syntaxin-1A", color=INK, fontsize=8.0, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, xQ + 0.26, Y + H - 0.25, 3)

# ---------------- Munc18-1 ---------------------------------------------------
mu = P(8.55, Y, 1.55, H)
mu.ladder("l", MU, {0, 1, 2, 3}, -1)
mu.draw(ax, SPEC, SPEC_ED, z=2)
ax.text(9.40, yM + 0.16, "Munc18-1", color=INK, fontsize=8.2, ha="center",
        va="center", fontweight="bold", zorder=6)
ax.text(9.40, yM - 0.26, "4/4", color=SPEC_ED, fontsize=7.6, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, 8.80, Y + H - 0.25, 3)

# ---------------- bridges ----------------------------------------------------
cx = P(3.35, Y + H, 2.20, 0.92).bottom(cpxV, -1, KNOB).bottom(cpxS, -1, KNOB)
cx.draw(ax, SPEC, SPEC_ED, z=2)
ax.text(4.55, Y + H + 0.64, "Complexin", color=INK, fontsize=7.8, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, 3.58, Y + H + 0.66, 2)

sy = P(5.55, Y + H, 2.15, 0.92).bottom(sytS, -1, KNOB).bottom(sytQ, -1, KNOB)
sy.draw(ax, SPEC, SPEC_ED, z=2)
ax.text(6.75, Y + H + 0.64, "Syt-1", color=INK, fontsize=7.8, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, 5.78, Y + H + 0.66, 2)

# ---------------- Synaptophysin ----------------------------------------------
sp = P(2.75, Y - 1.05, 2.35, 1.00)
sp.ladder("t", SP, SYP3, +1)
sp.draw(ax, SPEC, SPEC_ED, z=2)
ax.text(3.95, Y - 0.42, "Synaptophysin", color=INK, fontsize=7.6, ha="center",
        va="center", fontweight="bold", zorder=6)
ax.text(3.95, Y - 0.76, "3/5", color=SPEC_ED, fontsize=7.4, ha="center",
        va="center", fontweight="bold", zorder=6)
rings(ax, 3.00, Y - 0.28, 2)

# ---------------- SNCA -------------------------------------------------------
sn = P(1.05, Y, 1.95, H)
sn.ladder("r", NT, SYN2, +1)
sn.draw(ax, SYN, SYN_ED)
ax.text(1.90, yM + 0.18, "SNCA", color="#3a2c00", fontsize=9.2, ha="center",
        fontweight="bold", zorder=6)
ax.text(1.90, yM - 0.26, "2/5", color="#5c4700", fontsize=8.0, ha="center",
        fontweight="bold", zorder=6)
rings(ax, 1.30, Y + H - 0.25, 1, col="#5c4700")

# ---------------- socket markers ---------------------------------------------
ax.add_patch(plt.Circle((xV, yM), 0.33, facecolor="none", edgecolor=SYN_ED,
                        lw=1.6, ls=(0, (3, 2)), zorder=8))
ax.add_patch(plt.Circle((xS, yM), 0.33, facecolor="none", edgecolor=RETR_ED,
                        lw=1.6, ls=(0, (3, 2)), zorder=8))

# ---------------- ALTERNATIVE OCCUPANCY (not "losers") ------------------------
# AP180/CALM are not outcompeted: they bind the SAME socket in the OTHER state of
# the cycle, after NSF/alpha-SNAP disassembles the cis-SNARE complex. Overlap tells
# us they cannot bind SIMULTANEOUSLY; it does NOT tell us why. Competition and
# temporal succession give identical geometry.
ax.add_patch(FancyBboxPatch((5.35, 1.20), 4.45, 1.60, boxstyle="round,pad=0.06",
                            facecolor="white", edgecolor=RETR_ED, lw=1.3,
                            ls=(0, (5, 3)), zorder=1))
for i, nm in enumerate(["AP180", "CALM"]):
    x0 = 5.60 + i * 2.15
    pr = P(x0, 1.42, 1.85, 0.90)
    pr.ladder("l", [1.42 + 0.90 * f for f in (0.20, 0.35, 0.50, 0.65, 0.80)],
              ANTH3, +1)
    pr.draw(ax, RETR, RETR_ED)
    ax.text(x0 + 1.12, 1.98, nm, color=INK, fontsize=8.0, ha="center",
            fontweight="bold", zorder=6)
    ax.text(x0 + 1.12, 1.64, "3/5", color=STOP_ED, fontsize=7.6, ha="center",
            fontweight="bold", zorder=6)
    rings(ax, x0 + 0.25, 2.06, 3)
ax.add_patch(FancyArrowPatch((6.30, 2.90), (5.45, 4.30),
                             connectionstyle="arc3,rad=-0.24", arrowstyle="-|>",
                             mutation_scale=13, lw=1.4, color=RETR_ED,
                             ls=(0, (4, 2)), zorder=7))
ax.text(6.15, 3.62, "\u2297", color=STOP_ED, fontsize=12, ha="center",
        va="center", zorder=8)
ax.text(7.75, 3.12, "the same socket, the other state", color=RETR_ED,
        fontsize=6.9, ha="center", va="center", fontweight="bold", zorder=8)

fig.savefig("/home/claude/fig2_parkinsons.png", dpi=300, facecolor="white",
            bbox_inches="tight")
fig.savefig("/home/claude/fig2_parkinsons.pdf", facecolor="white", bbox_inches="tight")
print("wrote fig2 v2")
