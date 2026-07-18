from setuptools import setup, find_packages

with open("README.md", "r") as f:
    readme = f.read()

setup(
    name="gwdali",
    version="1.0.1",
    author="Josiel Mendonça Soares de Souza",
    author_email="josiel.souza@ufes.br",
    description="Upgrade of GWDALI with automatic-differentiation",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="BSD 3-Clause License",
    url="https://github.com/jmsdsouzaPhD/GWDALI/",
    keywords="fisher matrix gravitational waves gw dali jax",

    packages=find_packages(),
    include_package_data=True,

    package_data={
        "GWDALI": [
            "Sensitivities/*.txt",
            "*.png",
            "lib/*.txt",
            "lib/*.csv",
        ],
    },

    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "bilby",
        "astropy",
        "jax",
        "h5py",
    ],

    python_requires=">=3.9",

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)