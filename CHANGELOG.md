# Changelog

## v1.0.0

First public release; accompanies the manuscript submitted to *Bioinformatics Advances*.

**Two tiers**
- `scijigsaw-extract` — structures (PDB / mmCIF / predicted complexes) → interfaces, sites,
  overlaps. Contact cutoff, overlap threshold and confidence floor are user-set flags, not
  hidden defaults.
- `scijigsaw-render` — interaction table → laid-out, labelled jigsaw board (SVG/PDF/PNG).
  Figure 2 of the paper is produced by this command from the deposited CSVs.

**Assembly space**
- `scijigsaw-count` — exact linear-extension counting under the encoded seating model.
- `scijigsaw-bench` — 900-poset benchmark; regression and collinearity diagnostics.

**Reproducibility**
- `scripts/reproduce_numbers.py` regenerates every number in the paper and **exits non-zero
  if the code and the manuscript disagree**.
- `scripts/reproduce_figures.py` regenerates Figures 1–4 and S1.
- `scripts/benchmark_runtime.py` (optional, slow) measures the O(2ⁿ) runtime table.
- 24 tests; the paper's headline counts are pinned in `tests/test_assembly.py`.

**Known limitations, stated in the code as well as the paper**
- The seating model requires all partners present before a piece is seated. A multivalent
  protein may bind one partner and only later meet a second, so counts are a *lower bound*
  on biologically accessible pathways.
- Catalysis cannot be encoded: NSF/α-SNAP vacate the VAMP2 SNARE motif by a process with no
  static tab-and-socket representation.
- Contact coverage (n/N) proxies buried interface extent, **not** affinity.
- Bundled structures are **synthetic**, built to verify the extractor recovers known
  interfaces. Benchmarking against curated annotations on real complexes is **not yet done**.

## v1.0.1

Fixes three failures that only real structures expose. Each was invisible to the
synthetic test files, and each is now covered by `tests/test_structure_formats.py`.

- **A crystallographic B-factor is not a confidence score.** It is a temperature factor
  (~28 Å² for a well-ordered structure). Reading it as pLDDT refused every experimental
  structure. Predicted models are now detected from the file header; experimental
  structures are never refused for low B. New flag: `--confidence auto|plddt|none`.
- **Modified residues (MSE, SEP, TPO, PTR, …) are HETATM but are amino acids.** Skipping
  all HETATM deleted them from interfaces: six selenomethionines vanished from a
  sixteen-residue interface.
- **Insertion codes are distinct residues.** Residues are now keyed on (number, icode),
  so Kabat-numbered antibodies no longer conflate 52, 52A and 52B.

Also: `extract()` now reports `source` (predicted/experimental) alongside `confidence`.

## v1.1.0

**Parser robustness on deposited files.** Three failures that constructed complexes cannot
expose, each now covered by a regression test:

- a crystallographic **B-factor is a temperature factor, not a confidence score**. Read as
  pLDDT it falls below any sane floor, and the extractor refused every well-ordered
  experimental structure. Predicted models are now detected from the file header
  (`--confidence auto|plddt|none`).
- **modified residues** (MSE, SEP, TPO, PTR, …) are HETATM but are amino acids; discarding
  all HETATM deleted them from interfaces.
- **insertion codes** denote distinct residues (Kabat 52, 52A, 52B); keying on the residue
  number alone conflated them.

`examples/structures_parsing/` holds four real PDB entries used to test the parser. **They
are arbitrary structures, unrelated to the paper's biology, and validate nothing scientific.**

**Not done, and this is the substantive gap:** `scripts/validate_board_from_structures.py`
would apply the extractor to the deposited SNARE structures (1SFC, 3C98, 1KIL, 5CCG) and ask
whether the derived interfaces reproduce the board of Figure 2 — making it a derived result
rather than a curated one. It has not been run.

## v1.2.0 — the renderer can now draw the paper's own figure

The manuscript claimed Figure 2 was rendered from the deposited tables. **It was not.**
The renderer produced a hub-and-spoke board; the manuscript figure was a separate
hard-coded script. Since multivalent **bridge pieces** are the mechanic that does the
combinatorial work, a renderer that cannot draw a bridge cannot express the argument.

- **Schema extended** to general partner–partner edges, so SNAP25–Syntaxin and
  Syntaxin–Munc18 can be stated at all.
- **Topology is now derived, not declared.** The renderer finds the backbone (the chain
  of proteins binding ≥2 others), the bridges, the pendants and the contenders from the
  table alone.
- **A shared site is not always competition.** If two partners sharing a site on X are
  *adjacent to each other*, X lies in the groove between them — a bridge, not a contested
  site. Complexin's central helix contacts both VAMP2 and SNAP25 and is a bridge; AP180
  and SNAP25 both take VAMP2's SNARE motif and do **not** bind each other, so they
  compete. Getting this backwards benched half the board.
- Figure 2 of the paper is now produced by `scijigsaw-render` from `interactions.csv`.

41 tests, including `tests/test_topology.py`, which pins the derived backbone, bridges,
pendants and contenders.

## v1.3.0 — the 30S result is no longer hostage to a weak encoding

A referee's most damaging available line of attack was that the paper's headline contrast
rested on an encoding chosen to flatter it. The 30S hierarchy was reported under a
conservative "any-of" rule (a secondary requires *some* primary), giving **6×**. The
published Nomura/Held map names **specific** parents, and gives **45,545×** — a 7,500-fold
difference.

The conclusion survives, and is now stated to survive:

| assembly | subunits | depth | reduction |
|---|---:|---:|---:|
| NLRP3 inflammasome | 10 | 9 | 1,814,400× |
| 30S, *specific* parents (published map) | 20 | 3 | 45,545× |
| 30S, conservative "any-of" | 20 | 3 | 6× |

**Half the subunits, forty times the pruning**, even under the stronger encoding. The 30S
poset is only three deep either way. Depth, not size.

- `RIBOSOME_30S_SPECIFIC` added; both encodings are reported in the paper and in S4.3.
- `test_30S_result_is_not_an_artefact_of_a_weak_encoding` pins this: if a future change
  makes the contrast depend on the conservative choice, CI fails.

42 tests.
