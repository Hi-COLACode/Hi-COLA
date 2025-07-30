# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../../HiCOLA'))

import vinfo

# -- Project information

project = 'Hi-COLA'
# This needs to be changes to the Hi-COLA team.
# copyright = '2020-2025, Krishna Naidoo'
# author = 'Krishna Naidoo'

version = vinfo.vstr

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx_simplepdf'
]

autodoc_member_order = 'bysource'

source_suffix = ['.rst', '.md']

# Napoleon settings
napoleon_numpy_docstring = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

autodoc_typehints = 'none'

# -- Options for HTML output

#html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'

html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "HiCOLA_logo.jpg",
    "dark_logo": "HiCOLA_logo.jpg",
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
