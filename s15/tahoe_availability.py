#!/usr/bin/env python3
"""
Step 2 of the mTORC1 / rapamycin example.

Tahoe-100M supplies TWO things and nothing else:
  (i)  the perturbation condition, which indexes a declared model variant; and
  (ii) component availability, i.e. which subunits are expressed in a given
       cell line, which fixes the component SET of the board.

It supplies no contacts, no prerequisites and no exclusions. Expression is
never read as evidence about assembly order.

Asymmetry that governs the reporting: ABSENCE of transcript is strong evidence
(no transcript -> no protein -> component unavailable), whereas PRESENCE is weak
(transcript may not yield assembly-competent protein). Conclusions of the form
"this intermediate is IMPOSSIBLE in this context" are therefore robust;
conclusions of the form "formable" are permissive only. Report the negatives.

Usage
-----
    python tahoe_availability.py --out contexts.csv --threshold-sweep
"""
import argparse
import csv
import itertools
import sys

# mTORC1 subunit -> HGNC symbol. Verify against the Tahoe gene index; drop any
# subunit whose gene is not measured rather than silently treating it as absent.
SUBUNIT_GENES = {
    "mTOR":   "MTOR",
    "RAPTOR": "RPTOR",
    "mLST8":  "MLST8",
    "PRAS40": "AKT1S1",
    "DEPTOR": "DEPTOR",
    "FKBP12": "FKBP1A",
}

# Expression thresholds (pseudobulk CPM) to sweep. The reported default should be
# chosen the same way the 5 A contact cutoff was: state it, then show the
# conclusion is stable across the sweep (cf. Section S5.2).
THRESHOLDS = [0.5, 1.0, 2.0, 5.0]

DRUGS = ["rapamycin", "DMSO"]   # confirm the exact Tahoe compound label first


def load_pseudobulk(path):
    """Return {(cell_line, drug): {gene: cpm}}.

    Pseudobulk -- sum counts within each (cell_line, drug) group, then CPM
    normalise. Do NOT use per-cell values: single-cell dropout would manufacture
    false absences, and a false absence produces a false 'impossible' verdict,
    which is the one direction this analysis must not get wrong.
    """
    raise NotImplementedError(
        "Point this at the Tahoe-100M pseudobulk export. Keep dose, cell line "
        "and exposure time as provenance columns on every row."
    )


def component_set(expr, threshold):
    return {s for s, g in SUBUNIT_GENES.items()
            if expr.get(g, 0.0) >= threshold}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--out", default="contexts.csv")
    ap.add_argument("--threshold-sweep", action="store_true")
    args = ap.parse_args()

    data = load_pseudobulk(args.pseudobulk)
    thresholds = THRESHOLDS if args.threshold_sweep else [1.0]

    rows = []
    for (line, drug), expr in sorted(data.items()):
        if drug not in DRUGS:
            continue
        for t in thresholds:
            present = component_set(expr, t)
            rows.append({
                "cell_line": line,
                "drug": drug,
                "threshold_cpm": t,
                "model": "M_Rapa" if drug == "rapamycin" else "M0",
                "components_present": "|".join(sorted(present)),
                "components_absent": "|".join(
                    sorted(set(SUBUNIT_GENES) - present)),
            })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n_ctx = len({(r["cell_line"], r["drug"]) for r in rows})
    print(f"wrote {len(rows)} rows covering {n_ctx} (cell line, drug) contexts "
          f"-> {args.out}")
    print("Feed this to mtorc1_enumerate.py. Report only contexts where a "
          "subunit is ABSENT, and lead with the resulting impossibility calls.")


if __name__ == "__main__":
    main()
