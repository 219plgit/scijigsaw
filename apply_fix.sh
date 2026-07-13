#!/bin/bash
# Run this from INSIDE your local scijigsaw folder.
set -e

# 1. README: test count was stale (18 -> 24)
sed -i '' 's|pytest -q                    # 18 tests|pytest -q                    # 24 tests|' README.md

# 2. README: cite the release, not a DOI we are not minting
python3 - <<'PY'
import pathlib
r = pathlib.Path("README.md"); s = r.read_text()
s = s.replace("""## Citation

See `CITATION.cff`. Please cite both the software (archived DOI) and the paper.""",
"""## Citation

Cite the tagged release and the paper:

> Lio, P. and Lio, M.T. (2026) *scijigsaw: interface geometry as a constraint on
> protein-assembly order*, v1.0.0.
> https://github.com/219plgit/scijigsaw/releases/tag/v1.0.0

See `CITATION.cff`. A Zenodo DOI can be minted later by enabling the Zenodo-GitHub
hook and cutting a new release; nothing in the code needs to change.""")
r.write_text(s)

c = pathlib.Path("CITATION.cff"); t = c.read_text()
t = t.replace('doi: "INSERT-ZENODO-DOI"\n', 'url: "https://github.com/219plgit/scijigsaw/releases/tag/v1.0.0"\n')
t = t.replace('    orcid: "https://orcid.org/INSERT-ORCID"\n', '')
c.write_text(t)

p = pathlib.Path("pyproject.toml"); u = p.read_text()
u = u.replace('Archive    = "https://doi.org/INSERT-ZENODO-DOI"', 'Release    = "https://github.com/219plgit/scijigsaw/releases/tag/v1.0.0"')
p.write_text(u)
print("files updated")
PY

# 3. push scaffolding does not belong in a published tool
git rm -q --ignore-unmatch PUSH_INSTRUCTIONS.md

# 4. commit, MOVE the tag onto this commit, push both
git add -A
git commit -m "fix test count; cite the tagged release; drop push scaffolding"
git tag -f v1.0.0
git push origin main
git push -f origin v1.0.0

echo
echo "Done. Now publish the release:"
echo "  https://github.com/219plgit/scijigsaw/releases/new?tag=v1.0.0"
