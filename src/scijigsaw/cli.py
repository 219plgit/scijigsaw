"""Command-line entry points."""
from __future__ import annotations

import argparse
import json
import sys


def _assembly(argv=None):
    """scijigsaw-count : exact assembly-order counts for the encoded cases."""
    ap = argparse.ArgumentParser(
        description="Count assembly orders permitted by the encoded seating model.")
    ap.add_argument("case", choices=["vamp2", "inflammasome", "30S", "all"])
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)

    from .cases import ALL_CASES, count_30S
    out = {}
    if a.case in ("vamp2", "all"):
        out["vamp2"] = ALL_CASES["vamp2"].summary()
    if a.case in ("inflammasome", "all"):
        out["inflammasome"] = ALL_CASES["inflammasome"].summary()
    if a.case in ("30S", "all"):
        out["30S"] = count_30S()

    if a.json:
        print(json.dumps(out, indent=2))
        return
    for name, s in out.items():
        print(f"\n{name}")
        print(f"  units      {s['n']}")
        print(f"  depth      {s['depth']}")
        print(f"  orders n!  {s['total']:,}")
        print(f"  permitted  {s['permitted']:,}")
        print(f"  reduction  {s['reduction']:,.0f}x")
    print("\n  Counts are of orders permitted under the encoded SEATING model.")
    print("  They are a lower bound on biologically accessible pathways, not a")
    print("  characterisation of them (see assembly.py).")


def _bench(argv=None):
    """scijigsaw-bench : the random-poset benchmark."""
    ap = argparse.ArgumentParser(description="Random-poset benchmark.")
    ap.add_argument("--n-posets", type=int, default=900)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    from .benchmark import run, analyse
    A = run(a.n_posets, a.seed)
    r = analyse(A)
    print(f"posets: {r['n_posets']}\n")
    print("  Pearson r with log10(reduction)")
    for k, v in r["pearson"].items():
        print(f"    {k:<12} {v:+.3f}")
    print(f"\n  R2 (n, depth, width, density) = {r['r2_full']:.3f}")
    print(f"  R2 without depth              = {r['r2_without_depth']:.3f}")
    print(f"  depth adds                    = {r['depth_adds']:+.3f}")
    print("\n  variance-inflation factors")
    for k, v in r["vif"].items():
        print(f"    {k:<12} {v:5.1f}")
    print("\n  Depth, width and density are collinear: they are alternative")
    print("  measures of one property, the sequentiality of the assembly.")


def _extract(argv=None):
    """scijigsaw-extract : structures -> interaction table."""
    ap = argparse.ArgumentParser(description="Derive interfaces from structures.")
    ap.add_argument("structures")
    ap.add_argument("--out", default="interactions.csv")
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--overlap", type=float, default=0.20)
    ap.add_argument("--min-confidence", type=float, default=70.0)
    ap.add_argument("--min-interface-residues", type=int, default=3)
    ap.add_argument("--chains", default="A,B")
    ap.add_argument("--confidence", choices=["auto", "plddt", "none"], default="auto",
                    help="how to read the B-factor column. 'plddt': predicted model. "
                         "'none': experimental structure (B is a temperature factor, "
                         "NOT a confidence). 'auto': detect from the file header.")
    a = ap.parse_args(argv)

    from .extract import extract
    df, pairs = extract(a.structures, a.contact_cutoff, a.overlap,
                    a.min_confidence, a.min_interface_residues, a.chains)
    df.to_csv(a.out, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {a.out}")
    print(f"\n  contact cutoff {a.contact_cutoff} A, overlap {a.overlap}, "
          f"confidence floor {a.min_confidence}")
    print("  These are CHOICES. Interface detection is sensitive to them; sweep")
    print("  the cutoff and report the value used.")


def _render(argv=None):
    """scijigsaw-render : interaction table -> jigsaw board."""
    ap = argparse.ArgumentParser(description="Render a board from an interaction table.")
    ap.add_argument("proteins")
    ap.add_argument("interactions")
    ap.add_argument("--out", default="board.svg")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args(argv)

    from .render import render
    board, path = render(a.proteins, a.interactions, a.out, a.dpi)
    print(board.report())
    print(f"\nwrote {path}")


def _tiles(argv=None):
    """scijigsaw-tiles : interaction table -> printable cut-out kit(s)."""
    ap = argparse.ArgumentParser(
        description="Print an easy-to-cut set of tiles for building the complex by hand.")
    ap.add_argument("proteins")
    ap.add_argument("interactions")
    ap.add_argument("--out", default="tiles.pdf",
                    help="output file; .pdf gives multi-page A4 with an instruction/key "
                         "page, .svg a single stacked sheet")
    ap.add_argument("--variant", choices=["student", "teacher", "both"], default="both",
                    help="'student' = names only (a puzzle); 'teacher' = the answer key "
                         "with connector numbers, coverage and precedence/exclusion cues; "
                         "'both' writes _student and _teacher files")
    ap.add_argument("--title", default=None)
    a = ap.parse_args(argv)

    from .tiles import tiles
    kit, paths = tiles(a.proteins, a.interactions, a.out, a.variant, a.title)
    print(f"pieces     {len(kit.tiles)}")
    print(f"connectors {len(kit.connectors)}  "
          f"(each interface has its own connector shape; a tab fits only its socket)")
    for p in paths:
        print(f"wrote {p}")
    print("\n  Print at 100% and cut on the solid lines. The student set carries only")
    print("  protein names; matching connector shapes make it self-correcting. The")
    print("  teacher set adds the numbers, coverage and the precedence/exclusion answer.")
