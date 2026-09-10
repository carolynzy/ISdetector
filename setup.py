"""
File: setup.py
Project: ISdetector
Version: 1.1.0
Author: Yang Zhou
Purpose: Installation configuration for the ISdetector package.
"""
from setuptools import setup, find_packages
from pathlib import Path
from src import __version__

setup(
    name="ISdetector",
    version=__version__,
    description="Detect Insertion Sequences (IS) and associated Structural Variations (SVs) using paired-end or single-end WGS data.",
    long_description=Path(__file__).with_name("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Yang Zhou",
    packages=find_packages(),  # Automatically finds the 'src' package
    include_package_data=True,
    python_requires=">=3.10",
    # External executables cannot be enforced through Python package dependencies.
    install_requires=[
        "pysam>=0.19.0",
        "pandas",
        "biopython>=1.80",
        "numpy"
    ],
    entry_points={
        "console_scripts": [
            # Allows running the pipeline with the command 'isdetector'
            "isdetector=src.main:main",
        ],
    },
)
