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
