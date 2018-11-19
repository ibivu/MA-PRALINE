#!/usr/bin/env python

from __future__ import division, absolute_import, print_function

from setuptools import setup, find_packages, Extension
import numpy as np


entry_points = []
entry_points.append("PrositePatternAnnotator = mapraline.component:PrositePatternAnnotator")

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="ma-praline",
    version = "1.1.0",
    description = "PRALINE workflow for motif-aware alignment (MA-PRALINE).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author = "Maurits Dijkstra",
    author_email = "mauritsdijkstra@gmail.com",
    url="https://github.com/ibivu/MA-PRALINE/",
    license="GPLv2",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6"
    ],

    python_requires=">=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.2.*,!=3.3.*,!=3.4.*",
    install_requires = ["praline-aln>=1.1.0", "pyparsing>=2.0", "regex>=2018.11.07",
                        "six>=1.10.0"],
    packages = ["mapraline", "mapraline.component"],
    entry_points = {
      "praline.type": entry_points,
      "console_scripts": [
          "mapraline = mapraline.cmd:main",
          "mapraline_scores = mapraline.scores:main"
      ]
    }
)
