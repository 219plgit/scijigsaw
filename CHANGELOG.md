# Changelog

## v1.6.0

Structural-meaning and enumeration release accompanying the revised manuscript.
All reported counts are unchanged.

- **Enumerator.** The subset recurrence is now evaluated lazily, with maximality
  pruning, so only states reachable under the precedence and exclusion rules (the
  order ideals of the poset) are allocated. The VAMP2 board visits 40 of 128
  states; random posets of depth n/2 at n=20 visit 174 of 1,048,576. The worst
  case is unchanged at O(n*2^n) time and O(2^n) memory -- a poset with few
  relations still reaches every subset -- so the n=25 guard remains.
- **Renderer: n/N is no longer derived from affinity.** Previously the displayed
  n/N was computed from `kd_nM`, which contradicted the documented meaning. It is
  now taken only from explicit structural columns (`coverage`, `coverage_denom`)
  and is left blank for curated relations, so the literature-curated VAMP2 board
  displays no n/N. Core-clique selection now uses functional homogeneity and the
  alternative-occupancy occupant is chosen by interface degree; neither uses
  affinity.
- **Extractor.** New `--min-interface-residues` filter (default 3), matching the
  documented minimum interface size. New output fields
  `hub_interface_residue_count`, `overlap_group_union_residue_count`,
  `structural_footprint_fraction`, `coverage` and `coverage_denom`, so the
  structural footprint fraction |I_j|/|S_C| is computed from the same residue
  sets as the overlap coefficient.
- **Tiles and legend.** "binding sites" restated as "encoded interfaces";
  connector geometry described as symbolic ("a key, not a molecular surface");
  the n/N row restated structurally; fixed a right-edge text-clipping bug.
- **Figures.** Accessory connectors now seat flush with the core, and the core
  carries the tab throughout (bridges and accessories receive sockets), in both
  the board renderer and the assembly-steps script. Assembly-step labels resized
  and nudged to remove name overlaps; the in-figure payoff line moved to the
  caption; the benchmark figure overlays the three encoded assemblies.
- **Documentation.** The overlap coefficient |I_A n I_B| / min(|I_A|,|I_B|) and
  the single-linkage partner clustering are documented to match the code.

## v1.5.0

Physical output: build the complex with scissors.

- New `scijigsaw-tiles` prints an **easy-to-cut set of tiles** (multi-page A4 PDF,
  or a single SVG) in which every protein is a separate tile. Each interface is a
  complementary **tab and socket with its own keyed shape** (round / wedge /
  square, four sizes), so a tab fits only its true partner — the kit is
  self-correcting, a *control of error* in cardboard.
- Two variants from the same cut geometry: a **student** set (protein name only —
  a puzzle) and a **teacher** answer key that adds the connector number, coverage
  n/N, and the precedence (bridge, *seat last*) and exclusion (dashed,
  *either/or*) cues that `scijigsaw-count` enumerates, with an assembly-key page.
- Topology, colour and evolutionary rings are reused from the renderer; only the
  layout differs. Tests include a page-collision check (overlapping pieces cannot
  be cut) and a check that the student set carries no solution labels.

Reproducibility and robustness release accompanying the manuscript.

- Every value reported in the paper is pinned in the test suite and regenerated
  by `scripts/reproduce_numbers.py`, which exits non-zero on any mismatch; every
  figure by `scripts/reproduce_figures.py`.
- Random-poset benchmark (`scijigsaw-bench`): 900 posets with collinearity
  diagnostics (variance-inflation factors), showing that sequentiality — not unit
  count — governs how much of the order space geometry eliminates.
- Structure-extractor robustness on real deposited files (selenomethionine,
  alternate locations, insertion codes, multi-model files), with parsing tests.
- 51 tests across Python 3.11–3.13 under continuous integration.

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

## v1.5.0

- `scijigsaw.results.all_results()` — a single source of truth for every number the paper
  reports. The 336→252 correction had to be applied by hand in four documents and was missed
  in three; the manuscript, supplement, README and cover letter are now generated from one
  computed record.
- Test pinning the 1KIL structural correction, so it cannot silently regress.
