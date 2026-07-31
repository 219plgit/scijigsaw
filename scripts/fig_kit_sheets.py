#!/usr/bin/env python
"""Figures S3 and S4 -- the student and teacher tile sheets for one board.

Both sheets are drawn from the same TileKit with identical cut geometry, so a
teacher's solved model and a student's pieces interlock. The teacher variant
additionally carries connector identifiers, the bridge precedence cue and the
dashed alternative-occupancy outline.

    python scripts/fig_kit_sheets.py OUTDIR [proteins.csv interactions.csv]
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import pandas as pd

import scijigsaw.tiles as T
from scijigsaw.render import Board
from scijigsaw.tiles import TileKit

SUBTITLE = {
    "student": "All ten tiles on one sheet; the full printable kit is in the "
               "software repository.",
    "teacher": "The answer key: identical cut geometry to the student set, with "
               "connector numbers, precedence cues and dashed "
               "alternative-occupancy outlines.",
}


def sheet(kit, variant, out, board_name):
    tiles = [kit.tiles[n] for n in kit._pieces() if n in kit.tiles]
    A4_W, A4_H, M = T.A4_W, T.A4_H, T.MARGIN
    cols = 3
    rows = math.ceil(len(tiles) / cols)
    step_x = (A4_W - 2 * M) / cols
    step_y = T.TILE + 2 * T.CLEAR
    top = A4_H - M - 34
    pos = [(M + c * step_x + (step_x - T.TILE) / 2,
            top - (r + 1) * step_y + (step_y - T.TILE) / 2)
           for r in range(rows) for c in range(cols)]

    fig = plt.figure(figsize=(A4_W / 25.4, A4_H / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, A4_W); ax.set_ylim(0, A4_H)
    ax.set_aspect("equal"); ax.axis("off")
    ax.text(M, A4_H - M - 3, f"{board_name} board \u2014 {variant} tile set",
            fontsize=18, fontweight="bold", va="top")
    ax.text(M, A4_H - M - 15, SUBTITLE[variant], fontsize=11,
            color="#555555", va="top")
    for (x, y), t in zip(pos, tiles):
        kit._draw_tile(ax, t, x, y, variant)
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    prot = sys.argv[2] if len(sys.argv) > 2 else "examples/vamp2/proteins.csv"
    inter = sys.argv[3] if len(sys.argv) > 3 else "examples/vamp2/interactions.csv"
    os.makedirs(out, exist_ok=True)

    T.TILE, T.CLEAR = 33.0, 11.0            # sheet-figure scale
    board = Board(pd.read_csv(prot), pd.read_csv(inter))
    kit = TileKit(board)
    name = kit._hub_name() or "Encoded"
    sheet(kit, "student", os.path.join(out, "FigureS3_kit_student.pdf"), name)
    sheet(kit, "teacher", os.path.join(out, "FigureS4_kit_teacher.pdf"), name)


if __name__ == "__main__":
    main()
