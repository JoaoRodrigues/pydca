# Copyright 2016 by Joao Rodrigues.  All rights reserved.
# This code is part of the pydca distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Test suite for HMMERWrapper.
"""

import os
import sys
import unittest

from pydca.io import HMMERWrapper

class WrapperTest(unittest.TestCase):

    def setUp(self):
        """Common across all tests"""

        self.hw = HMMERWrapper.HMMERWrapper

        modpath = os.path.abspath(os.path.dirname(__file__))
        self.seqfile = os.path.join(modpath, 'data', 'P00929.fasta')
        self.badfile = os.path.join(modpath, 'data', 'bad.fasta')

    def test_fromfile(self):
        """Reading FASTA file"""

        wrapper = self.hw(self.seqfile, mock=True)
        n_seqs = len(wrapper.sequences)
        self.assertEqual(n_seqs, 2)

    def test_fromhandle(self):
        """Reading file-like object"""

        with open(self.seqfile, 'r') as handle:
            wrapper = self.hw(handle, mock=True)

    def test_fromlist(self):
        """Reading an unsupported input format (and failing)"""

        self.assertRaises(TypeError, self.hw, [])

    def test_readbadformat(self):
        """Reading a bad FASTA file (and failing)"""

        self.assertRaises(HMMERWrapper.ParseError, self.hw, self.badfile)

if __name__ == '__main__':
    unittest.main(verbosity=2)
