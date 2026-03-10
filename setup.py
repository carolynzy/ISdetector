"""
File: setup.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Installation configuration for the ISdetector package.
"""
from setuptools import setup, find_packages
setup(
    name="ISdetector",
    version="1.0.0",
    description="Detect Insertion Sequences (IS) and associated Structural Variations (SVs) using paired-end or single-end WGS data.",
    author="Yang Zhou",
    packages=find_packages(),  # Automatically finds the 'src' package
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "pysam>=0.19.0",   # 
        "pandas",          # 
        "biopython>=1.80", # 
        "numpy"
    ],
    entry_points={
        "console_scripts": [
            # Allows running the pipeline with the command 'isdetector'
            "isdetector=src.main:main",
        ],
    },
)
