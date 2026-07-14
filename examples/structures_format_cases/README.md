# Structural-format edge cases

**These files are SYNTHETIC.** They are not experimental structures and they are not
real PDB entries. They are hand-built to exercise features that every real structure has
and that the toy complexes in `../structures/` do not — nothing more.

| file | what it exercises | the bug it caught |
|---|---|---|
| `EXPT__partner.pdb` | `CRYST1` + `EXPDTA`; B-factors ≈ 28 Å² | a crystallographic B-factor is a **temperature factor, not a confidence score**. Read as pLDDT it falls below any sane floor, and the extractor refused *every* experimental structure. |
| `SEMET__partner.pdb` | selenomethionine (MSE) as `HETATM` at the interface | modified residues are HETATM but **are** amino acids. Skipping all HETATM deleted six MSE residues from a sixteen-residue interface. |
| `ICODE__partner.pdb` | insertion codes (52, 52A, 52B, 52C) | Kabat-numbered antibodies carry insertion codes. Keying on the residue number alone **conflates distinct residues**. |

Passing these tests demonstrates that the parser handles these formats. It does **not**
demonstrate accuracy on experimental structures: interface-residue precision and recall
against curated annotations remains unmeasured. See `../../README.md`, *Status of
validation*.
