#!/usr/bin/env python3
"""
Immunoproteasome subunit switching as a Scientific Jigsaw model index.

WHY THIS EXAMPLE IS CENTRAL RATHER THAN PERIPHERAL
--------------------------------------------------
The earlier proteasome example used perturbation as an external WEIGHT. Under
Plackett-Luce that can never change the admissible set:

    Omega_control = Omega_drug,    only p(omega) moves

Here the perturbation changes COMPONENT AVAILABILITY, which changes the
component set, which changes Omega itself:

    Omega_control != Omega_drug

Three catalytic beta-ring positions each admit exactly one of two mutually
exclusive occupants (a declared EXCLUSION, the third relation type and the one
least exercised elsewhere in the manuscript):

    position     constitutive      inducible
    beta1        PSMB6             PSMB9  (beta1i, LMP2)
    beta2        PSMB7             PSMB10 (beta2i, MECL-1)
    beta5        PSMB5             PSMB8  (beta5i, LMP7)

Evidence for the exclusions is STRUCTURAL (one subunit per ring position),
never transcriptomic. Tahoe supplies only which variants are available.

Note the systematic numbering does not track the gene numbering:
beta1 = PSMB6, beta5 = PSMB5, beta6 = PSMB1. State this in the manuscript.

"Intermediate proteasomes" carrying mixed constitutive and inducible subunits
are documented in normal tissue and in human cancer cell lines, so which of the
eight compositions are formable is a real open question, not a hypothetical --
and it is exactly a Jigsaw support query.

This board is HUMAN by necessity (the immunoproteasome is vertebrate-specific),
which removes the yeast/human transfer assumption that burdened the S14.2
board.

    python3 immunoproteasome_board.py
    python3 immunoproteasome_board.py --availability tahoe_availability.csv
"""
import argparse
import csv
import itertools
import os
import sys

# beta-ring position -> (constitutive gene, inducible gene)
POSITIONS = {
    "b1": ("PSMB6", "PSMB9"),
    "b2": ("PSMB7", "PSMB10"),
    "b5": ("PSMB5", "PSMB8"),
}
CONSTITUTIVE = {c for c, _ in POSITIONS.values()}
INDUCIBLE = {i for _, i in POSITIONS.values()}

# Non-catalytic beta subunits, present in every composition. Included so the
# board is a real 20S half-ring rather than only the catalytic triple; they do
# not vary with the perturbation and so do not affect the composition count.
STRUCTURAL = ["PSMB1", "PSMB2", "PSMB3", "PSMB4"]


def compositions(available):
    """Every composition consistent with the declared exclusions.

    Exactly one occupant per position. A position with neither variant
    available cannot be filled, so the complex is not formable at all --
    the enumerator returns nothing, which is an exact impossibility verdict,
    not a low score.
    """
    choices = []
    for pos, (c, i) in POSITIONS.items():
        opts = [v for v in (c, i) if v in available]
        if not opts:
            return [], pos          # unfillable position named
        choices.append(opts)
    return list(itertools.product(*choices)), None


def orders_per_composition(n_catalytic=3):
    """Admissible join orders for one composition.

    PLACEHOLDER: 3! = 6, i.e. the three catalytic subunits join a 13S-like
    seed in any order. REPLACE with the project enumerator once the human
    contact graph is extracted (PDB 5LE5 / 4R3O or equivalent), so that
    declared contacts and prerequisites constrain this properly.
    """
    import math
    return math.factorial(n_catalytic)


def analyse(label, available):
    comps, blocked = compositions(available)
    if blocked:
        return {"condition": label, "compositions": 0, "omega": 0,
                "verdict": f"IMPOSSIBLE: position {blocked} unfillable",
                "mixed": 0}
    per = orders_per_composition()
    mixed = sum(1 for c in comps
                if set(c) & CONSTITUTIVE and set(c) & INDUCIBLE)
    return {"condition": label, "compositions": len(comps),
            "omega": len(comps) * per,
            "verdict": "formable", "mixed": mixed}


def scenarios():
    C, I = CONSTITUTIVE, INDUCIBLE
    return [
        ("control: constitutive only", C),
        ("full switch: inducible only", I),
        ("mixed: all six available", C | I),
        ("partial: b5i induced", C | {"PSMB8"}),
        ("partial: b5i + b1i induced", C | {"PSMB8", "PSMB9"}),
        ("PSMB6 lost, no b1i available", (C | I) - {"PSMB6", "PSMB9"}),
        ("PSMB5 lost, b5i available", (C | I) - {"PSMB5", "PSMB9", "PSMB10"}),
    ]


def observed_queries(available):
    """Support queries: is a specific reported composition admissible?

    Intermediate proteasomes reported in the literature go here. Each is a
    query against the declared model, never an input to it.
    """
    named = {
        "constitutive (b1,b2,b5)":      ("PSMB6", "PSMB7", "PSMB5"),
        "immuno (b1i,b2i,b5i)":         ("PSMB9", "PSMB10", "PSMB8"),
        "intermediate b1/b2/b5i":       ("PSMB6", "PSMB7", "PSMB8"),
        "intermediate b1/b2i/b5i":      ("PSMB6", "PSMB10", "PSMB8"),
    }
    comps, blocked = compositions(available)
    cset = {frozenset(c) for c in comps}
    return {k: ("admissible" if frozenset(v) in cset else "not admissible")
            for k, v in named.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--availability",
                    help="CSV with columns cell_line,drug,gene,expressed")
    ap.add_argument("--out", default="immunoproteasome_spaces.csv")
    a = ap.parse_args()

    print("Board: human 20S catalytic beta ring")
    print(f"  positions   : {', '.join(f'{p}({c}|{i})' for p,(c,i) in POSITIONS.items())}")
    print(f"  structural  : {', '.join(STRUCTURAL)} (invariant)")
    print(f"  exclusions  : {len(POSITIONS)} declared, one occupant per position")
    print(f"  orders/composition: {orders_per_composition()}  "
          f"[PLACEHOLDER -- replace with extracted contact graph]\n")

    rows = []
    if a.availability:
        import collections
        avail = collections.defaultdict(set)
        with open(a.availability) as fh:
            for r in csv.DictReader(fh):
                if str(r["expressed"]).strip().lower() in ("1", "true", "yes"):
                    avail[(r["cell_line"], r["drug"])].add(r["gene"])
        for (line, drug), genes in sorted(avail.items()):
            rows.append(analyse(f"{line} / {drug}", genes))
    else:
        rows = [analyse(l, s) for l, s in scenarios()]

    w = max(len(r["condition"]) for r in rows) + 2
    print(f"{'condition':<{w}}{'comps':>7}{'|Omega|':>9}{'mixed':>7}  verdict")
    for r in rows:
        print(f"{r['condition']:<{w}}{r['compositions']:>7}{r['omega']:>9}"
              f"{r['mixed']:>7}  {r['verdict']}")

    if not a.availability:
        print("\nSupport queries under 'mixed: all six available':")
        for k, v in observed_queries(CONSTITUTIVE | INDUCIBLE).items():
            print(f"  {k:<28} {v}")
        print("\nUnder 'control: constitutive only':")
        for k, v in observed_queries(CONSTITUTIVE).items():
            print(f"  {k:<28} {v}")

    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nwrote {a.out}")
    print("\nThe control -> mixed change is a change of Omega, not of p(omega).")
    print("That is what distinguishes model INDEXING from external WEIGHTING.")
    os._exit(0)


if __name__ == "__main__":
    main()
