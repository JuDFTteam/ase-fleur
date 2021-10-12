# -*- coding: utf-8 -*-
# Copyright (c) Materials Virtual Lab
# Distributed under the terms of the Modified BSD License.
"""
setup: usage: pip install -e .
"""
import os
from setuptools import setup, find_packages


SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(SETUP_PTH, "README.md")) as f:
    desc = f.read()


setup(
    name="ase-fleur",
    packages=find_packages(),
    version="0.0.1",
    install_requires=[
        "ase",
        "masci-tools>=0.5.0",
    ],
    entry_points={
        "ase.ioformats": ["fleur-inpgen = ase_fleur.io:inpgen_format", "fleur-xml= ase_fleur.io:xml_format"],
        "ase.calculator": ["fleur = ase_fleur.calculator:Fleur"],
    },
    extras_require={},
    package_data={},
    author="Henning Janssen",
    author_email="he.janssen@fz-juelich.de",
    maintainer="Henning Janssen",
    url="https://github.com/janssenhenning/ase-fleur",
    license="BSD",
    description="Package adding IO/calculator functionalities for the FLEUR code to the ase package.",
    long_description=desc,
    keywords=["ase"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
