# Deposited structures — used as PARSER tests only

These are real PDB entries, but they are **not** structures of anything in the paper.
They were chosen because they exercise features of deposited files that constructed
complexes cannot: crystallographic B-factors, cryo-EM and NMR provenance, non-`A`/`B`
chain identifiers, multiple models, HETATM records.

| PDB | method | what it exercises |
|---|---|---|
| 2XHE | X-ray | B-factors are temperature factors, not confidence scores |
| 7DDO | cryo-EM | the partner chain is `C`, not `B` |
| 2BEG | NMR | Aβ(1–42) fibril; a control — a fully buried layer must return every residue |
| 1LCD | NMR | multiple models; protein–DNA |

**They validate the parser. They validate nothing about the science.** They are not
structures of VAMP2, SNAP25, syntaxin, Munc18 or complexin, and no scientific claim in
the paper rests on them.

The test that *would* validate the science is `scripts/validate_board_from_structures.py`,
which asks whether the extractor, applied to the deposited structures of the SNARE
assemblies themselves (1SFC, 3C98, 1KIL…), recovers the board that Figure 2 asserts.
**That test has not been run.**
