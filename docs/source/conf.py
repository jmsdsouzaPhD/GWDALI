# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'GWDALI'
copyright = '2024, Josiel'
author = 'Josiel Mendonça Soares de Souza'

release = '0.0'
version = '0.1.4'

master_doc = 'index'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme_options = {
    'style_nav_header_background': 'linear-gradient(180deg, rgba(200,200,200,1) 50%, rgba(0,0,0,1) 100%)',
}
html_theme = 'sphinx_rtd_theme'
#html_theme = 'groundwork'
html_logo = 'logo_gwdali_name.png'

# -- Options for EPUB output
epub_show_urls = 'footnote'
