#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup
from glob import glob


__version__ = "0.0.1-dev"


classes = """
    Development Status :: 2 - Pre-Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 3.4
    Programming Language :: Python :: 3.5
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

long_description = ("Horizomer: workflow for whole genome HGT detection.")

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='horizomer',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Horizomer',
      author="Horizomer development team",
      author_email="jenya.kopylov@gmail.com",
      url='https://github.com/biocore/horizomer',
      test_suite='nose.collector',
      packages=['benchmark', 'distance-method'],
      package_data={},
      scripts=glob('benchmark/*py') + glob('distance-method/*py') +
      glob('benchmark/tests/*py'),
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'doc': ["Sphinx == 1.3.3"]},
      install_requires=['click >= 6'],
      classifiers=classifiers
      )
