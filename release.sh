#!/usr/bin/env bash
# scijigsaw v3.0.0 — Phase 1 (GitHub release for submission)
# Run from the repository root AFTER merging the final modules and
# applying the metadata files from this package. Requires: git, gh (authenticated).
set -euo pipefail

VERSION="3.0.0"
TAG="v${VERSION}"

# 0. sanity: clean tree, tests, version agreement
git diff --quiet || { echo "uncommitted changes — commit first"; exit 1; }
python -m pytest -q
grep -q "version = \"${VERSION}\"" pyproject.toml || { echo "pyproject version mismatch"; exit 1; }
grep -q "version: ${VERSION}" CITATION.cff || { echo "CITATION.cff version mismatch"; exit 1; }

# 1. reproduce the manuscript numbers (should be green before tagging)
python reproduce_all.py

# 2. set today's date in CITATION.cff
sed -i "s/^date-released: .*/date-released: $(date +%F)/" CITATION.cff
git add CITATION.cff
git commit -m "v${VERSION}: set release date" || true

# 3. annotated tag + push
git tag -a "${TAG}" -m "scijigsaw ${TAG}: evidence-weighted inference and hypothesis navigation (PLOS submission release)"
git push origin main --follow-tags

# 4. GitHub release from the prepared notes
gh release create "${TAG}" \
  --title "scijigsaw ${TAG} — evidence-weighted inference and hypothesis navigation" \
  --notes-file RELEASE_NOTES_v3.0.0.md

# 5. optional: attach the submission archive as an asset
# gh release upload "${TAG}" scijigsaw_${TAG}_submission_archive.tar.gz

echo "Done. Release: https://github.com/219plgit/scijigsaw/releases/tag/${TAG}"
echo "Phase 2 (Zenodo + DOI) happens after acceptance — see RELEASE_CHECKLIST."
