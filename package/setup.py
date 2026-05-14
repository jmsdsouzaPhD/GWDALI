from setuptools import setup, find_packages

with open("../README.md",'r') as arq:
	readme = arq.read()

setup(
	name = 'gwdali',
	version = '1.0',
	license = 'BSD 3-Clause License',
	author  = 'Josiel Mendonça Soares de Souza',
	author_email = 'josielsouza@if.ufrj.br',
	keywords = 'fisher matrix, gravitational waves, gw, dali, jax',
	description = 'Upgrade of GWDALI with automatic-differentiation',
	packages = find_packages(),
	include_package_data=True,
	package_data={
        'GWDALI': ['Sensitivities/*.txt'],
    },
	install_requires = [
						'numpy',
						'matplotlib',
						'scipy',
						'bilby',
						'astropy',
						#'jax'
						],
	url = "https://github.com/jmsdsouzaPhD/GWDALI/",
)
