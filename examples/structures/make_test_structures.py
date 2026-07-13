"""Generate the synthetic two-chain complexes used to TEST the extractor.

These are NOT real structures. They are constructed with KNOWN interfaces so we
can verify that the extractor recovers what it is given -- a correctness check,
not evidence of accuracy on real data. Benchmarking against curated interface
annotations on real complexes remains to be done (see the paper, Section 10.1).

Ground truth built in:
    SNAP25        hub residues 30-50   (the SNARE motif)
    AP180         hub residues 30-40   -> OVERLAPS SNAP25, must be flagged
    SNCA          hub residues  1-12   at 6.5 A -> absent at 5 A, present at 8 A
    Synaptophysin hub residues 55-60
    Complexin     hub residues 32-44 but partner pLDDT 30 -> must be REFUSED
"""
import numpy as np


def atom(i, name, res, chain, rnum, xyz, b):
    return ("ATOM  %5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n"
            % (i, name, " ", res, chain, rnum, xyz[0], xyz[1], xyz[2], b))


def write(fn, span, plddt, gap=4.0):
    lines, i = [], 1
    for r in range(1, 61):                       # chain A = hub, 60 residues
        lines.append(atom(i, "CA", "ALA", "A", r, (r * 3.8, 0.0, 0.0), 90.0)); i += 1
    for k, r in enumerate(span):                 # chain B = partner
        lines.append(atom(i, "CA", "GLY", "B", k + 1, (r * 3.8, gap, 0.0), plddt)); i += 1
    open(fn, "w").writelines(lines + ["END\n"])


if __name__ == "__main__":
    write("VAMP2__SNAP25.pdb",        range(30, 51), 95)
    write("VAMP2__AP180.pdb",         range(30, 41), 92)
    write("VAMP2__SNCA.pdb",          range(1, 13),  88, gap=6.5)
    write("VAMP2__Synaptophysin.pdb", range(55, 61), 85)
    write("VAMP2__Complexin.pdb",     range(32, 45), 30)
    print("wrote 5 synthetic complexes with known interfaces")
