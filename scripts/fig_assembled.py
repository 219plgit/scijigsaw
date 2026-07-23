"""Clean assembled-complex schematic (the teacher answer-key scheme, Figure S4).
Single clean tab/socket per interface; functional colours and names only.
"""
import sys, os, pandas as pd, matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType (no Type 3)
matplotlib.rcParams['ps.fonttype'] = 42; import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from collections import defaultdict
from scijigsaw.render import Board, INK
from scijigsaw.geometry import Piece, TAB, SOCKET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else 'FigureS4_assembled.pdf'
b = Board(pd.read_csv(os.path.join(ROOT,'examples/vamp2/proteins.csv')),
          pd.read_csv(os.path.join(ROOT,'examples/vamp2/interactions.csv')))
bb = b.backbone
W, H, y0, x0 = 2.0, 1.5, 4.2, 0.8
xs = {p: x0 + i*W for i,p in enumerate(bb)}
yM = y0 + H/2
R = 0.24                                   # absolute connector radius (world units)
kn = lambda edge_len: R/edge_len           # -> knob as a fraction of that edge

def style(n):
    fc, ec = b._style(n); return fc, ec

# ---- plan bridge feet (match core top tabs) ----
bridge_feet = {}; core_top = defaultdict(list)
for br,(a,c) in b.bridges.items():
    if bb.index(a) > bb.index(c): a,c = c,a
    xa, xc = xs[a] + W*0.70, xs[c] + W*0.30
    bridge_feet[br] = (xa, xc)
    core_top[a].append(xa); core_top[c].append(xc)

# ---- plan pendants (match core bottom sockets) ----
by_host = defaultdict(list)
for pd_, host in b.pendants.items(): by_host[host].append(pd_)
pend = {}; core_bot = defaultdict(list); GAP=0.18
for host, pds in by_host.items():
    pds = sorted(pds); k=len(pds); span=W*0.84; pw=(span-(k-1)*GAP)/k; sx=xs[host]+(W-span)/2
    for j,pd_ in enumerate(pds):
        px = sx + j*(pw+GAP); tabX = px + pw/2
        pend[pd_] = (px, pw, tabX); core_bot[host].append(tabX)

fig, ax = plt.subplots(figsize=(9.2, 6.4), facecolor='white')
ax.set_aspect('equal'); ax.axis('off')
ax.set_xlim(0.2, x0 + len(bb)*W + 0.3)
ax.set_ylim(y0 - 1.18 - 0.35, y0 + H + 0.92 + 0.45)

# ---- core row ----
for i,p in enumerate(bb):
    pc = Piece(xs[p], y0, W, H)
    if i>0: pc.left(yM, SOCKET, knob=kn(H))
    if i<len(bb)-1: pc.right(yM, TAB, knob=kn(H))
    for X in core_top[p]: pc.top(X, TAB, knob=kn(W))
    for X in core_bot[p]: pc.bottom(X, TAB, knob=kn(W))
    fc,ec = style(p)
    ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.8, zorder=3))
    ax.text(xs[p]+W/2, yM, p, color=INK, fontsize=10.5, ha='center', va='center', fontweight='bold', zorder=6)

# ---- bridges above ----
for br,(xa,xc) in bridge_feet.items():
    bx0, bx1 = min(xa,xc)-0.36, max(xa,xc)+0.36
    bw = bx1-bx0; by=y0+H; bh=0.92
    pc = Piece(bx0, by, bw, bh)
    pc.bottom(xa, SOCKET, knob=kn(bw)); pc.bottom(xc, SOCKET, knob=kn(bw))
    fc,ec = style(br)
    ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.7, zorder=2))
    ax.text((bx0+bx1)/2, by+bh/2+0.16, br, color=INK, fontsize=9.5, ha='center', va='center', fontweight='bold', zorder=6)

# ---- accessories below ----
ph = 1.18
for pd_, (px,pw,tabX) in pend.items():
    pc = Piece(px, y0-ph, pw, ph)
    pc.top(tabX, SOCKET, knob=kn(pw))
    fc,ec = style(pd_)
    ax.add_patch(Polygon(pc.polygon(), closed=True, facecolor=fc, edgecolor=ec, lw=1.7, zorder=2))
    ax.text(px+pw/2, y0-ph/2-0.05, pd_, color=INK, fontsize=8.8, ha='center', va='center', fontweight='bold', zorder=6)

fig.savefig(OUT, facecolor='white', bbox_inches='tight')
print("clean assembled board written")
