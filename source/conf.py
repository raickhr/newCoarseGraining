# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = 'coarseGraining'
copyright = '2025, Shikhar Rai'
author = 'Shikhar Rai'
release = '0.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

# Support both reStructuredText and Markdown files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
