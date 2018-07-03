#!/usr/bin/env python

from __future__ import division, absolute_import, print_function

from setuptools import setup, find_packages, Extension
import numpy as np


entry_points = []
entry_points.append('PrositePatternAnnotator = mapraline.component:PrositePatternAnnotator')

setup(name='Motif-aware PRALINE',
      version = '1.1',
      description = 'PRALINE workflow for motif-aware alignment (MA-PRALINE).',
      author = 'Maurits Dijkstra',
      author_email = 'm.j.j.dijkstra@vu.nl',
      url = 'http://www.few.vu.nl/',
      license = "GPL",

      install_requires = ['PRALINE>=1.1', 'pyparsing>=2.0', 'regex>=2.4.46', 'six>=1.10.0'],
      packages = ['mapraline', 'mapraline.component'],
      entry_points = {
          'praline.type': entry_points,
          'console_scripts': [
              'mapraline = mapraline.cmd:main',
              'mapraline_scores = mapraline.scores:main'
          ]
      }

    )
