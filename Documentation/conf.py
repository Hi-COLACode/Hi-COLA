"""Build guides and API reference from the checked-out Hi-COLA source."""
import ast
import os
from pathlib import Path
import sys

ROOT = Path.cwd() if (Path.cwd() / "HiCOLA").is_dir() else Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
project = "Hi-COLA"
author = "Hi-COLA team"
release = os.environ.get("DOCS_VERSION", "latest")
version = release
extensions = [
    "myst_parser", "sphinx.ext.autodoc", "sphinx.ext.autosummary",
    "sphinx.ext.napoleon", "sphinx.ext.mathjax",
]
root_doc = "index"
exclude_patterns = ["_build", "_generated", "README.md", "FML-dependency-audit.md"]
autosummary_generate = True
autodoc_typehints = "description"
napoleon_numpy_docstring = True
myst_heading_anchors = 3
myst_enable_extensions = ["dollarmath", "colon_fence"]
html_theme = "sphinx_book_theme"
templates_path = ["_templates"]
repository = os.environ.get("GITHUB_REPOSITORY", "Hi-COLACode/Hi-COLA")
source_ref = os.environ.get("DOCS_SOURCE_REF", "Hi-COLA_public")
html_theme_options = {
    "logo": {
        "image_light": "_static/Hi-COLA-light.png",
        "image_dark": "_static/Hi-COLA-dark.png",
        "text": "Hi-COLA",
    },
    "repository_url": f"https://github.com/{repository}",
    "repository_branch": source_ref,
    "path_to_docs": "Documentation",
    "use_repository_button": True,
    "use_source_button": True,
    "home_page_in_toc": True,
}
# Left sidebar: logo, version selector, search, then the page list.
html_sidebars = {
    "**": ["navbar-logo", "release-switcher", "search-button-field", "sbt-sidebar-nav"],
}
# The repository logos are copied into _static alongside the theme assets.
html_static_path = ["_static", "../logos"]
html_favicon = "../logos/Hi-COLA-light.png"
html_js_files = ["release-switcher.js"]
myst_url_schemes = {
    "http": None, "https": None, "mailto": None,
    "repo": {"url": f"https://github.com/{repository}/blob/{source_ref}/{{{{path}}}}"},
}

# Enumerate definitions without importing executable CLI modules.
modules = {
    "Expression builder": "HiCOLA.Frontend.expression_builder",
    "Numerical solver": "HiCOLA.Frontend.numerical_solver",
    "Parameter reader": "HiCOLA.Frontend.read_parameters",
    "Utilities": "HiCOLA.Utilities.Other.support",
    "Backend launcher": "HiCOLA.backend",
    "Plotting: power spectrum": "HiCOLA.Plotting.pofk",
    "Plotting: background": "HiCOLA.Plotting.background",
    "Plotting: growth": "HiCOLA.Plotting.growth",
    "Plotting: frontend tables": "HiCOLA.Plotting.frontend",
    "Plotting: shared helpers": "HiCOLA.Plotting.common",
}
lines = ["Python API", "==========", "",
         "Signatures and docstrings are generated from the checked-out source.",
         "Functions without docstrings show their signatures only.", ""]
for title, module in modules.items():
    source = ROOT / (module.replace(".", "/") + ".py")
    # Modules absent from this checkout are omitted from the API reference.
    if not source.is_file():
        continue
    tree = ast.parse(source.read_text())
    names = [node.name for node in tree.body
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
             and not node.name.startswith("_")]
    lines.extend([title, "-" * len(title), "", ".. autosummary::",
                  "   :toctree: generated", "   :nosignatures:", ""])
    lines.extend(f"   {module}.{name}" for name in names)
    lines.append("")
(ROOT / "Documentation/api").mkdir(parents=True, exist_ok=True)
(ROOT / "Documentation/api/index.rst").write_text("\n".join(lines))

# Example figures are generated for each build; see figures.py.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures  # noqa: E402

figures.generate(ROOT, ROOT / "Documentation")
