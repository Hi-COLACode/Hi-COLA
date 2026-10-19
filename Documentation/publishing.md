# Building and publishing documentation

[Documentation home](index.md)

The site uses Sphinx, MyST Markdown, and the PyData theme. Python API pages come from source definitions and docstrings. Executable generator and plotting modules are excluded from imports. Use NumPy-style docstrings for new functions.

## Local build

Use Python 3.12.14 and the pinned documentation dependencies:

```sh
python -m pip install -r requirements-docs.lock
MPLBACKEND=Agg python -m sphinx -W --keep-going -b html Documentation Documentation/_build/html
lychee --offline --include-fragments --no-ignore --exclude-path '(^|/)_static/' -- 'Documentation/_build/html/**/*.html'
```

Generated API stubs and build output are ignored by Git. The dropdown reads the shared `versions.json` file on the published site; a standalone local HTML build does not include that manifest.

## Publication and GitHub Pages

`.github/workflows/documentation.yml` validates documentation on pull requests. Relevant default-branch pushes and manual runs build only the current documentation, then commit the assembled site to **`gh_pages`** using GitHub Pages Deploy Action. Existing version directories are retained, while `current/` is replaced. Git history is preserved.

The workflow reads GitHub's default-branch setting, currently `Hi-COLA_public`, so renaming that branch needs no workflow edit. Documentation, frontend/API source, dependencies, packaging, README/logo assets, and workflow changes trigger builds. Superseded runs for the same branch or PR are cancelled; publication jobs are serialised to avoid concurrent site updates. Pull requests never publish.

The workflow explicitly deploys the committed `gh_pages` content using `actions/upload-pages-artifact` and `actions/deploy-pages`. A branch push using `GITHUB_TOKEN` does not itself trigger GitHub's automatic branch-based Pages build. Select **GitHub Actions** in repository **Settings → Pages**. Allow default-branch and release-tag deployments in the `github-pages` environment. Publication requires `contents: write`, `pages: write`, and `id-token: write`; no additional secret is needed.

The site root redirects to `current/`. A shared version list drives the dropdown, including on archived pages. Switching versions keeps the current page when it exists and otherwise opens the selected version's home page.

## Release archives

Publishing a GitHub release triggers archiving under its tag name. If the tag matches the saved current build, `current/` is moved into the release directory before current documentation is replaced. If the newly built current docs match, they are copied into the release directory. Otherwise, that tag alone is checked out and built once using its own lock file. Historical releases are not rebuilt on ordinary default-branch pushes.

Each directory includes `source-commit.txt`. Existing release directories are not overwritten. Republishing a release with an unchanged tag reuses the archive; moving an already archived tag fails publication. Release tags must be single directory components containing letters, digits, dots, underscores, or hyphens, and cannot use reserved names such as `current`. Tags requiring a build must contain the Sphinx setup and documentation lock file. GitHub releases, rather than every Git tag, create archives; existing historical tags are not automatically backfilled.

Current source links point to its exact source commit, so archived copies retain the correct source links without editing their HTML. The Python version and direct/transitive package versions are pinned. Pinning versions does not guarantee identical bytes across operating systems or toolchains; the generated release HTML is retained unchanged in `gh_pages`.

## Independent link workflow

`.github/workflows/links.yml`, named **Links**, uses Lychee for both local and external checks. Local checks run offline and validate files and anchors in rendered HTML. Theme templates are excluded from scan inputs; assets referenced by pages are still checked. Make **Links / local-links** a required branch-protection check to prevent merging broken local links.

External checks scan documentation, rendered HTML, and Python, C++, Lua, and INI source files directly, including URLs in comments, docstrings, strings, and code blocks. Failures are advisory because external sites may be unavailable or reject automated requests. Lychee adds a job summary and saves `documentation-link-report`; local results are saved as `local-link-report`.
