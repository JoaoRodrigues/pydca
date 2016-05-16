# Copyright 2016 by Joao Rodrigues.  All rights reserved.
# This code is part of the pydca distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Distutils based setup script for pydca

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard Python mechanism
for installing packages. For the easiest installation just type the command:

    python setup.py install

For more in-depth instructions, see the installation section of the manual.

If all else fails, feel free to write to the author(s) and ask for help.
"""

from __future__ import print_function

from distutils.core import setup
import sys

# Check the Python version
if sys.version_info[:2] < (2, 7):
    print('pydca requires Python 2.7 or later. '
          'Python {0.major}.{0.minor} detected'.format(sys.version_info))
    sys.exit(1)

setup(name='pydca',
      version='0.0.1',
      description='Coevolution analysis for protein sequences',
      author='Joao Rodrigues, Jason Wang, Jessica Zhao',
      author_email='j.p.g.l.m.rodrigues@gmail.com',
      url='http://github.com/drorlab/pydca',
      packages=['pydca', 'pydca.core', 'pydca.io', 'pydca.wrappers'],
     )
