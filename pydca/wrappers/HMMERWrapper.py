# Copyright 2016 by Joao Rodrigues.  All rights reserved.
# This code is part of the pydca distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Wrapper class to produce a single Stockholm file with custom annotations from a
jackhmmer/phmmer sequence search.

Depends heavily on Biopython SearchIO and AlignIO modules.
"""

from __future__ import print_function

from collections import namedtuple
import logging
import os
import tempfile

try:
    # Python 3.3+
    from shutil import which as _which
except ImportError:
    def _which(executable):
        """
        Returns None if the executable is not found in the PATH.
        Source: http://stackoverflow.com/a/377028
        """

        def is_exe(fpath):
            """Returns True if the path is an executable"""
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, _ = os.path.split(executable)
        if fpath:
            if is_exe(executable):
                return executable
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, executable)
                if is_exe(exe_file):
                    return exe_file

        return None

try:
    from Bio import SeqIO
    from Bio.Alphabet.IUPAC import IUPACProtein
    from Bio.Alphabet import _verify_alphabet
except ImportError, err:
    raise ImportError('[!] Could not import one or more biopython modules: {}'.format(err))

# Constants
CANONICAL_AA1 = set('ACDEFGHIKLMNPQRSTVWY')

# Classes
class HMMERRuntimeError(Exception):
    """Custom class to throw exceptions related to running HMMER"""
    pass

class ParseError(Exception):
    """Custom class to throw exceptions related to parsing input/output"""
    pass

class HMMERWrapper(object):
    """Wraps phmmer/jackhmmer to produce a custom-annotated Stockholm file.

    Uses Biopython to call HMMER (default jackhmmer) on one or more protein sequences (file
    or handle) and parses several output files to create a unified Stockholm file with the
    full aligned sequences and a set of features for each hit (e.g. uniprot accession, hit
    e-value, etc) and a set of per-file annotations (e.g. number of sequences, sequence
    length, original HMMER parameters, etc). This file is a suitable input to create an
    Alignment object.

    Args:
        sequence: sequence file, open file handle, or string with the data in FASTA format
        database: a string with the path to the sequence database against which to run HMMER

        executable: a string with the HMMER executable to call (optional) [def: jackhmmer]
        evalue: a number in scientific notation with the expectation value threshold (optional)
        ncpus: an integer with the number of CPUs to use in the HMMER calculation (optional)
        niter: an integer with the number of iterations to use in jackhmmer (optional)

        cleanup: a boolean to remove all files produced by HMMER
        mock: a boolean to block the actual call to HMMER (use for testing only)

    Returns:
        stockof: the path to the Stockholm output file.

    Raises:
        HMMERRuntimeError: An error occured while executing HMMER.
        ParseError: An error occured while parsing the output of HMMER.
        OSError: HMMER executable cannot be found.
        TypeError: input/output does not match expected format.
    """

    def __init__(self, sequence, database=None, executable='jackhmmer', **kwargs):

        # Empty container for parsed/validated sequences
        self.sequences = []

        self.database = database

        # Setup logging with module name
        self.logger = logging.getLogger(__name__)

        # Iterate over kwargs and define defaults if not user-provided
        _defaults = {'evalue': 1E-20,
                     'ncpus': 2,
                     'niter': 5,
                     'mock': False}

        for kwarg in _defaults:
            kwarg_val = kwargs.get(kwarg)
            if not kwarg_val:
                kwarg_val = _defaults[kwarg]
            setattr(self, kwarg, kwarg_val)

        # Validate sequence data and locate executable
        self.__validate(sequence)

        if not _which(executable):
            raise OSError("HMMER executable not found in PATH: '{}'".format(executable))


        # run HMMER

    def __validate(self, seqdata):
        """Verifies the user provided input is either a file or a file-like object
        containing sequence data.

        Args:
            seqdata: a string (path to file) or file-like object.

        Returns:
            A list with one namedtuple per input sequence.

        Raises:
            TypeError: if the input is not a string or a file-like object.
            ParseError: if the sequence contains others than the 20 canonical AAs.
        """

        _Sequence = namedtuple('Seq', ['name', 'data'])

        # file-like object
        # isinstance(obj, file) does not hold in Py3
        if hasattr(seqdata, 'read') and hasattr(seqdata, 'name'):
            self.logger.debug('Reading data from file-like object {}'.format(seqdata.name))
            fname = seqdata.name

        elif isinstance(seqdata, basestring):
            self.logger.debug('Reading data from file path {}'.format(seqdata))
            fname = seqdata

            # can be file name string or sequence
            if not os.path.isfile(fname):
                raise OSError('Sequence file not found: {}'.format(seqdata))
        else:
            raise TypeError('Sequence input format not recognized: {}'.format(seqdata))

        # parse and validate sequences
        # defining these two a prior just in case later we decide to support more stuff
        _seq_alphabet = IUPACProtein()
        _seq_format = 'fasta'

        seq_iterator = SeqIO.parse(seqdata, _seq_format, alphabet=_seq_alphabet)
        for seq_i, seq_record in enumerate(seq_iterator, start=1):

            seq_name = seq_record.name
            seq_raw = str(seq_record.seq)
            if not _verify_alphabet(seq_record.seq):
                msg = 'Entry #{} ({}) in {} is not a valid protein sequence'
                raise ParseError(msg.format(seq_i, seq_name, fname))

            self.sequences.append(_Sequence(seq_name, seq_raw))

        return self.sequences

    def run(self):
        """Launches the HMMER executable on the user data"""
        pass

