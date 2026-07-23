# VAMP2-centred real-structure benchmark

A compact, systematic evaluation of the structure-extraction tier on the deposited
complexes that underlie the VAMP2 board. It turns the single 1KIL correction into a
reproducible experiment with precision/recall, a cutoff sweep, and the resulting change
in the encoded constraint graph.

## Structures (download into `structures/`)

| PDB  | Context                                   | Tests                                              |
|------|-------------------------------------------|----------------------------------------------------|
| 1SFC | Core neuronal SNARE complex               | recovery of VAMP2–Syntaxin–SNAP25 core interfaces  |
| 1KIL | SNARE–complexin complex                   | correct complexin bridge; reject complexin–SNAP25  |
| 5W5C | Primed SNARE–complexin–synaptotagmin-1    | multivalent regulatory interfaces                  |
| 3C98 | Munc18a–syntaxin-1 complex                | Munc18-1–Syntaxin edge                             |

Optional robustness replicates: **5W5D** (same study as 5W5C — a replicate, not an
independent system) and **3RK2** (truncated prefusion SNARE–complexin).

Download, e.g.:

    mkdir -p structures
    for id in 1SFC 1KIL 5W5C 3C98; do
      curl -L "https://files.rcsb.org/download/${id}.cif" -o "structures/${id}.cif"
    done

## Level 1 — chain-pair interface detection (runnable here)

Reference positive/negative chain pairs are curated in `reference_labels.csv`
(**verify each against the cited structure paper before reporting**). Then:

    python run_benchmark.py --structures ./structures --labels reference_labels.csv \
        --cutoffs 4 5 6 8 --min-residues 3 5 10 --out results.csv

`run_benchmark.py` runs the real extractor at each parameter setting and reports
precision, recall and F1 per structure and in aggregate. Treat **5 Å / 5 residues** as
the predefined default; the others are the sensitivity analysis. Do not tune a separate
cutoff per structure.

## Level 2 — residue-level agreement (needs an external reference)

For chain pairs correctly called as interfaces, compare the residues returned by the
extractor with the interface residues reported by **PDBePISA** or the **RCSB pairwise
interface** record, and report the Jaccard index J = |A∩B| / |A∪B| and residue-level
precision/recall. State openly that this is *agreement with an established interface
resource*, not independent biological validation. The PISA/RCSB residue sets are not
bundled; fetch them per structure and place them in `reference_residues/`.

## Consequence for the encoded graph

For each parameter setting, also record the number of detected interfaces, the number of
false interfaces, and the resulting **VAMP2-board order count**. This links extractor
performance to the paper's central output (e.g. whether small residue changes leave the
precedence graph — and therefore the 252/5,040 count — unchanged, or whether crossing a
threshold changes it). Report whichever behaviour occurs; do not assume it.

## Scope

The benchmark assesses conversion of deposited structures into chain-pair, residue-level
and shared-site constraints. It does not validate assembly chronology, kinetic selection,
or literature-only relations (AP180/CALM– and α-synuclein–VAMP2) for which no suitable
complex structure is available.
