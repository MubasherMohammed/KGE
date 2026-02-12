"""
Setup script for MAGeCK-KGE (MAGeCK with interactive HTML reports).

This installs the modified mageck Python package on top of the
bioconda-installed mageck, preserving the compiled binaries (mageckGSEA)
while replacing the Python source with our enhanced version.
"""

from setuptools import setup, find_packages

setup(
    name='mageck-kge',
    version='0.5.9.5.kge1',
    description='MAGeCK with interactive Plotly HTML reports and pathway enrichment',
    author='Mubasher Mohammed',
    packages=find_packages(),
    package_data={
        'mageck': [
            '*.gmt',
            '*.txt',
            '*.Rnw',
            '*.Rmd',
            '*.RTemplate',
        ],
    },
    entry_points={
        'console_scripts': [
            'mageck=mageck.cli:main',
        ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'scikit-learn',
        'matplotlib',
        'plotly>=6.0',
    ],
    python_requires='>=3.11',
)
