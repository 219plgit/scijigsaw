# Changelog

## v1.4.0

Applied the extractor to the deposited complexin–SNARE complex (PDB 1KIL) and corrected
the encoding it contradicted.

- `Complexin` bridges **VAMP2 and Syntaxin-1A**, not VAMP2 and SNAP25; it makes no contact
  with SNAP25. Added the direct VAMP2–Syntaxin edge (34 interface residues): the SNARE core
  is a four-helix bundle, not a chain.
- **Permitted orders 336 → 252 of 5,040** (93.3% → 95.0% eliminated).
- Core detection rewritten: a bundle is a clique, so the core is the most tightly bound
  clique, with a longest-path fallback for nucleated cascades.

## v1.3.0

- The 30S hierarchy is reported under both encodings. Conservative "any-of": 6×. The
  specific parents named by the published map: 45,545×. The inflammasome (10 subunits,
  depth 9) still prunes 40× more than the 30S (20 subunits, depth 3) — depth, not size.
- `RIBOSOME_30S_SPECIFIC` added, with a test pinning the contrast.

## v1.2.0

- The renderer now draws the paper's Figure 2 from `interactions.csv`. It previously
  produced only a hub-and-spoke board.
- Schema extended to partner–partner edges; backbone, bridges, pendants and contenders are
  derived from the table rather than declared.
- A shared site means competition only when the two partners do **not** bind each other; if
  they do, the piece bridges them.

## v1.1.0

Parser fixes for deposited structures — none of which synthetic test files can expose.

- A crystallographic B-factor is a temperature factor, not a confidence score. Reading it as
  pLDDT refused every experimental structure. Predicted models are now detected from the
  header (`--confidence auto|plddt|none`).
- Modified residues (MSE, SEP, TPO, PTR…) are HETATM but are amino acids; discarding all
  HETATM deleted them from interfaces.
- Insertion codes denote distinct residues (Kabat 52, 52A, 52B).

## v1.0.0

First public release.

- `scijigsaw-extract` — structures → interfaces, sites, overlaps. Contact cutoff, overlap
  threshold and confidence floor are user-set flags.
- `scijigsaw-render` — interaction table → jigsaw board.
- `scijigsaw-count` — exact linear-extension counting under the encoded seating model.
- `scijigsaw-bench` — 900-poset benchmark.
- `scripts/reproduce_numbers.py` regenerates every number in the paper and exits non-zero if
  the code and the manuscript disagree; `scripts/reproduce_figures.py` regenerates the figures.

## Known limitations

- The seating model requires all partners present before a piece is seated. A multivalent
  protein may bind one partner and meet a second later, so counts are a lower bound on
  biologically accessible pathways.
- Catalysis cannot be encoded: NSF/α-SNAP vacate the VAMP2 SNARE motif by an ATP-dependent
  process with no static tab-and-socket representation.
- Contact coverage (n/N) proxies buried interface extent, not affinity.
- The extractor has not been benchmarked against curated interface annotations.
