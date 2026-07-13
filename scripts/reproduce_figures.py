#!/usr/bin/env python3
"""Regenerate every FIGURE in the paper and supplement.  ~40 s.

    python scripts/reproduce_figures.py --out figures/

Produces (PDF + PNG):
    Figure1_tool            the rule, four channels, two tiers
    Figure2_VAMP2           the VAMP2 board                 [from the CSVs, via the renderer]
    Figure3_depth           three encoded assemblies
    Figure4_benchmark       900 random posets
    FigureS1_inflammasome   the NLRP3 board: a chain, not a hub
"""
import argparse, math, os, subprocess, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scijigsaw import INFLAMMASOME, VAMP2, count_30S, render
from scijigsaw.benchmark import run

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INK, MUTED = "#1a1d24", "#5b6675"


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300,
                    facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def fig3_depth(out):
    rows = [("NLRP3 inflammasome", INFLAMMASOME.summary(), "#1f7a6f"),
            ("VAMP2 board", VAMP2.summary(), "#b8860b")]
    r30 = count_30S()
    fig, ax = plt.subplots(figsize=(6.6, 4.3), facecolor="white")
    for name, s, col in rows:
        ax.scatter(s["depth"], s["reduction"], s=130, color=col, zorder=3,
                   edgecolor="white", linewidth=1.4)
        ax.annotate(f"{name}\n(n={s['n']})", (s["depth"], s["reduction"]),
                    textcoords="offset points", xytext=(12, -4), fontsize=8, color=INK)
    ax.scatter(r30["depth"], r30["reduction"], s=130, color="#c0392b", zorder=3,
               edgecolor="white", linewidth=1.4)
    ax.annotate(f"30S ribosomal subunit\n(n={r30['n']})",
                (r30["depth"], r30["reduction"]), textcoords="offset points",
                xytext=(12, -4), fontsize=8, color=INK)
    ax.set_yscale("log"); ax.set_xlim(1.5, 11.5)
    ax.set_xlabel("depth of the dependency poset  (longest chain)", fontsize=10)
    ax.set_ylabel("assembly orders eliminated  (n! / permitted)", fontsize=10)
    ax.grid(alpha=.25, which="both", lw=.6)
    ax.tick_params(labelsize=9, colors=MUTED)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#c7cdd6")
    fig.tight_layout(); save(fig, out, "Figure3_depth")


def fig4_benchmark(out):
    A = run(n_posets=900, seed=0)
    n_, d_, w_, dens_, br_, y = A.T
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), facecolor="white")
    ax = axes[0]
    rng = np.random.default_rng(0)
    sc = ax.scatter(n_ + rng.uniform(-.18, .18, len(n_)), y, c=d_, s=16,
                    cmap="viridis", alpha=.85, edgecolor="none")
    cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label("poset depth", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    ax.set_xlabel("number of units, n", fontsize=10)
    ax.set_ylabel("orders eliminated,  log$_{10}$(n! / L)", fontsize=10)
    ax.set_title("(a)  size alone does not determine pruning", fontsize=9.5, color=INK)
    ax = axes[1]
    for n0, col in [(10, "#9ecae1"), (12, "#4292c6"), (14, "#08519c")]:
        m = n_ == n0
        ds = sorted(set(d_[m].astype(int)))
        ax.plot(ds, [y[m & (d_ == d)].mean() for d in ds], "o-", color=col,
                ms=4.5, lw=1.6, label=f"n = {n0}")
    ax.set_xlabel("poset depth (longest chain)", fontsize=10)
    ax.set_ylabel("mean orders eliminated,  log$_{10}$", fontsize=10)
    ax.set_title("(b)  depth does, at every n", fontsize=9.5, color=INK)
    ax.legend(fontsize=8.5, frameon=False)
    for ax in axes:
        ax.grid(alpha=.25, lw=.6); ax.tick_params(labelsize=8.5, colors=MUTED)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"): ax.spines[sp].set_color("#c7cdd6")
    fig.tight_layout(); save(fig, out, "Figure4_benchmark")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print("Figure 2 -- the VAMP2 board, rendered FROM THE CSVs")
    render(os.path.join(ROOT, "examples/vamp2/proteins.csv"),
           os.path.join(ROOT, "examples/vamp2/interactions.csv"),
           os.path.join(a.out, "Figure2_VAMP2.pdf"))
    render(os.path.join(ROOT, "examples/vamp2/proteins.csv"),
           os.path.join(ROOT, "examples/vamp2/interactions.csv"),
           os.path.join(a.out, "Figure2_VAMP2.png"))
    print("  wrote Figure2_VAMP2.pdf / .png")

    print("Figure 3 -- depth, not size")
    fig3_depth(a.out)
    print("Figure 4 -- the random-poset benchmark")
    fig4_benchmark(a.out)

    for script, name in [("fig1_min.py", "Figure1_tool"),
                         ("fig4_inflammasome.py", "FigureS1_inflammasome")]:
        p = os.path.join(HERE, script)
        if os.path.exists(p):
            print(f"{name} -- {script}")
            subprocess.run([sys.executable, p, os.path.abspath(a.out)],
                           check=True, capture_output=True)
            print(f"  wrote {name}.pdf / .png")
    print(f"\nAll figures written to {a.out}/")


if __name__ == "__main__":
    main()
