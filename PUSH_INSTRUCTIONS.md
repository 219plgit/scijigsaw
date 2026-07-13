# Publishing this repository

The git history, tag `v1.0.0` and remote are already set up. You need only push.

## 1. Create the empty repo on GitHub

Go to https://github.com/new
  - Owner: **219plgit**
  - Name:  **scijigsaw**
  - Public
  - **Do not** add a README, .gitignore or licence (they already exist here)

## 2. Push

```bash
cd scijigsaw
git push -u origin main
git push origin v1.0.0
```

If you use SSH rather than HTTPS:

```bash
git remote set-url origin git@github.com:219plgit/scijigsaw.git
```

## 3. Connect Zenodo (this is what mints the DOI)

1. Sign in at https://zenodo.org with your GitHub account.
2. Go to **Settings → GitHub**, find `219plgit/scijigsaw`, and switch it **On**.
3. Back on GitHub: **Releases → Draft a new release**
   - Tag: `v1.0.0` (it already exists — just select it)
   - Title: `scijigsaw v1.0.0`
   - Description: paste from `CHANGELOG.md`
   - **Publish release**
4. Zenodo archives it automatically and mints a DOI. `.zenodo.json` supplies the
   metadata (title, authors, affiliations, licence, keywords), so nothing needs
   editing there.

Note: Zenodo only archives releases created **after** the switch is turned on. If you
publish the release first and enable the hook afterwards, delete the release and
recreate it.

## 4. Paste the DOI in five places

Zenodo gives you two DOIs. Use the **concept DOI** (the "all versions" one) — it always
resolves to the latest release.

| file | placeholder |
|---|---|
| `pyproject.toml` | `Archive = "https://doi.org/INSERT-ZENODO-DOI"` |
| `CITATION.cff` | `doi: "INSERT-ZENODO-DOI"` |
| `README.md` | add a DOI badge (Zenodo gives you the markdown) |
| manuscript | *Data and software availability* section |
| manuscript | the software reference: `Lio,P. and Lio,M.T. (2026) scijigsaw ... doi:[INSERT-ZENODO-DOI]` |
| supplement | S9 |

Then commit the DOI and (optionally) cut `v1.0.1`.

## 5. Check CI is green

`.github/workflows/test.yml` runs the 24 tests on Python 3.9 / 3.11 / 3.12. The README
badge will go green once it passes. **If CI is red, do not submit** — the tests pin the
paper's numbers, so a red badge means the code and the manuscript disagree.

## Still outstanding for the paper (not for the repo)

The extractor has been validated only against **synthetic** structures with known
interfaces. Benchmarking it against curated interface annotations on 20–30 real complexes
— precision, recall, site-clustering agreement, cutoff sensitivity — remains the one
substantive scientific gap, and both the README and the manuscript say so.
