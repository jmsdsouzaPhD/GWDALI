import os
import sys

sys.path.insert(0, os.path.abspath('../../package'))

project = 'GWDALI'
copyright = '2026, Josiel'
author = 'Josiel Mendonça Soares de Souza'

release = '1.0'
version = '1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

autosummary_generate = True

templates_path = ['_templates']
html_static_path = ['_static']

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only': True,
    'style_nav_header_background': 'linear-gradient(180deg, rgba(200,200,200,1) 50%, rgba(0,0,0,1) 100%)',
}

html_logo = '_static/logo_gwdali.png'