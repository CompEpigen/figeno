import os
import sys
import mock
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

MOCK_MODULES = ['pyBigWig', 'pysam', 'numpy', 'pandas', 'vcfpy', 'scipy','statsmodels','pyyaml','matplotlib','cooler','Flask','importlib-resources']
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()

project = 'figeno'
copyright = '2024, Etienne Sollier'
author = 'Etienne Sollier'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx_rtd_theme','sphinx.ext.autosectionlabel','sphinxarg.ext']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
html_theme = 'sphinx_rtd_theme'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
