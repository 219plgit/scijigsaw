"""Teacher assembly-order exercise for the VAMP2 board (Figure S4).

Three progressive panels -- close the core seam, seat the bridges, add the
accessories -- that let a class build the complex by hand and discover the
precedence and alternative-occupancy constraints the software enumerates.
Single clean tab/socket per interface; functional colours; TrueType (no Type 3).
"""
import sys, os, pandas as pd, matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from collections import defaultdict
from scijigsaw.render import Board, INK
from scijigsaw.geometry import Piece, TAB, SOCKET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else 'FigureS4_assembly.pdf'
b = Board(pd.read_csv(os.path.join(ROOT, 'examples/vamp2/proteins.csv')),
          pd.read_csv(os.path.join(ROOT, 'examples/vamp2/interactions.csv')))
bb = b.backbone
W, H, y0, x0 = 2.0, 1.5, 4.2, 0.8
xs = {p: x0 + i * W for i, p in enumerate(bb)}
yM = y0 + H / 2
R = 0.24
kn = lambda e: R / e
style = lambda n: b._style(n)

bridge_feet = {}; core_top = defaultdict(list)
for br, (a, c) in b.bridges.items():
    if bb.index(a) > bb.index(c): a, c = c, a
    xa, xc = xs[a] + W * 0.70, xs[c] + W * 0.30
    bridge_feet[br] = (xa, xc); core_top[a].append(xa); core_top[c].append(xc)
by_host = defaultdict(list)
for pd_, host in b.pendants.items(): by_host[host].append(pd_)
pend = {}; core_bot = defaultdict(list); GAP = 0.18
for host, pds in by_host.items():
    pds = sorted(pds); k = len(pds); span = W * 0.84
    pw = (span - (k - 1) * GAP) / k; sx = xs[host] + (W - span) / 2
    for j, pd_ in enumerate(pds):
        px = sx + j * (pw + GAP); tabX = px + pw / 2
        pend[pd_] = (px, pw, tabX, host); core_bot[host].append(tabX)

XL = (0.2, x0 + len(bb) * W + 0.3)
YL = (y0 - 1.18 - 0.30, y0 + H + 0.92 + 0.30)

def draw_board(ax, level):
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(*XL); ax.set_ylim(*YL)
    for i, p in enumerate(bb):
        pc = Piece(xs[p], y0, W, H)
        if i > 0: pc.left(yM, SOCKET, knob=kn(H))
        if i < len(bb) - 1: pc.right(yM, TAB, knob=kn(H))
        if level >= 2:
            for X in core_top[p]: pc.top(X, TAB, knob=kn(W))
        if level >= 3:
            for X in core_bot[p]: pc.bottom(X, TAB, knob=kn(W))
        fc, ec = style(p)
        ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.8, zorder=3))
        ax.text(xs[p] + W / 2, yM, p, color=INK, fontsize=7.2, ha='center', va='center', fontweight='bold', zorder=6)
    if level >= 2:
        for br, (xa, xc) in bridge_feet.items():
            bx0, bx1 = min(xa, xc) - 0.36, max(xa, xc) + 0.36; bw = bx1 - bx0; by = y0 + H; bh = 0.92
            pc = Piece(bx0, by, bw, bh); pc.bottom(xa, SOCKET, knob=kn(bw)); pc.bottom(xc, SOCKET, knob=kn(bw))
            fc, ec = style(br)
            ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.7, zorder=2))
            ax.text((bx0 + bx1) / 2, by + bh / 2 + 0.16, br, color=INK, fontsize=7.2, ha='center', va='center', fontweight='bold', zorder=6)
    if level >= 3:
        ph = 1.18
        for pd_, (px, pw, tabX, host) in pend.items():
            pc = Piece(px, y0 - ph, pw, ph); pc.top(tabX, SOCKET, knob=kn(pw))
            fc, ec = style(pd_)
            ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.7, zorder=2))
            cx = px + pw / 2; hc = xs[host] + W / 2
            nudge = 0.08 if cx > hc + 1e-6 else (-0.08 if cx < hc - 1e-6 else 0.0)
            ax.text(cx + nudge, y0 - ph / 2 - 0.05, pd_, color=INK, fontsize=6.5, ha='center', va='center', fontweight='bold', zorder=6)

fig, axes = plt.subplots(3, 1, figsize=(6.6, 9.7), facecolor='white')
steps = [(1, "Step 1  \u2014  close the core seam",
             "SNAP25\u2013Syntaxin-1A\u2013VAMP2 join first. The seam must close before anything can seat on top."),
         (2, "Step 2  \u2014  seat the bridges",
             "Syt-1 and Complexin each need two sockets, so they attach only once the seam beneath them is closed (precedence)."),
         (3, "Step 3  \u2014  add the accessories",
             "Each accessory joins its own partner. On VAMP2's SNARE-motif site, SNAP25, AP180 or CALM bind \u2014 one at a time (alternative occupancy).")]
for ax, (lvl, title, desc) in zip(axes, steps):
    draw_board(ax, lvl)
    ax.text(XL[0], YL[1] + 0.05, title, fontsize=11.5, fontweight='bold', color=INK, va='bottom', ha='left')
    ax.text(XL[0], YL[0] - 0.02, desc, fontsize=8.0, color='#444444', va='top', ha='left', wrap=True)
fig.subplots_adjust(top=0.95, bottom=0.03, hspace=0.42)
fig.savefig(OUT, facecolor='white', bbox_inches='tight')
print("stepwise assembly exercise written:", OUT)
