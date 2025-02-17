# -*- coding: utf-8 -*-




import os
import sys
import roshambo2
import pkg_resources


sys.path.append(os.path.abspath("../"))


# version specified in ../setup.py
version = pkg_resources.require("roshambo2")[0].version



extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "myst_parser",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "member-order": "bysource",
}

autoclass_content = 'both'  # Can be 'both', 'class', or 'init'
autodoc_member_order = 'bysource'  # Preserves the order of members


source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

master_doc = "index"

project = "Roshambo2"
copyright = "2025"


exclude_patterns = ["_build", "_templates"]
# html_static_path = ["_static"]
# templates_path = ["_templates"]

pygments_style = "sphinx"

html_theme = "pydata_sphinx_theme"


# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
