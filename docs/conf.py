# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../"))


# -- Project information -----------------------------------------------------

project = "gufe"
copyright = "2022, The OpenFE Development Team"
author = "The OpenFE Development Team"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinxcontrib.autodoc_pydantic",
    "sphinx.ext.intersphinx",
]

autoclass_content = "both"

autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    # Need to stop at ints to avoid docs breaking:
    # https://stackoverflow.com/questions/77914856/
    # this originates from the docstrings of to_bytes
    # and from_bytes in cpython.
    "inherited-members": "GufeTokenizable, int",
    "undoc-members": True,
}

autosummary_generate = True

intersphinx_mapping = {
    "rdkit": ("https://www.rdkit.org/docs/", None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_mock_imports = [
    "msgpack",
    "rdkit",
]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
# html_theme_options = {
#    "accent_color": "FeelingFabulous",
# }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# replace macros
rst_prolog = """
.. |rdkit.mol| replace:: :class:`rdkit.Chem.rdchem.Mol`
"""
