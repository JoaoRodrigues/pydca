# Copyright 2016 by Joao Rodrigues, Jason Wang, Jessica Zhao.  All rights reserved.
# This code is part of the pydca distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Data structures for multiple sequence alignments.
"""

from __future__ import division

from collections import namedtuple
import logging

import numpy as np

# Classes
class AlignmentError(Exception):
    """Custom class to throw exceptions related to creating and modfying Alignment objects"""
    pass

# Lightweight Sequence object
Sequence = namedtuple('Sequence', ['name', 'data', 'annotations'])

class Alignment(object):
    """
    Class to represent multiple sequence alignment.
    Implements add/remove/lookup methods and generic alignment statistics methods.
    """

    def __init__(self, seq_data, seqnames, seq_annotations=None, aln_annotations=None):
        """
        Initializes an Alignment object.

        Args:
            seq_data [list]: list of lists (2D matrix) with one sequence per 'row'
            seqnames [list]: Unique identifiers for each sequence.
            seq_annotations [dict]: Dictionary with per-sequence/generic annotations.
                                    If key matches sequence name, added to seq annotation.
            aln_metadata [dict]: high-level alignment annotations (GC/GF entries).
        """

        self.logger = logging.getLogger(__name__)

        n_seq = len(seq_data)
        n_snm = len(seqnames)
        if n_seq != n_snm:
            emsg = 'Different numbers of sequences and sequence identifiers: {0} vs {1}'
            raise AlignmentError(emsg.format(n_seq, n_snm))

        self.sequences = np.array(seq_data, dtype='S1')
        self.logger.debug('Initialized 2D array of shape {0}'.format(self.sequences.shape))
        self.seq_annotations = {seq_name: {} for seq_name in seqnames}
        self.aln_annotations = aln_annotations

        self.annotate(seq_annotations)

        # sequence ID to index mapping and vice-versa
        self.names = seqnames
        self._mapping_r = {seq_n: seq_i for seq_i, seq_n in enumerate(self.names)}

    # Private Methods
    def __repr__(self):
        """String representation of the Alignment object"""
        _qu, (_len, _aa) = self.query, self.sequences.shape
        return "<Alignment of {0} with {1} sequences of length {2}>".format(_qu, _len, _aa)

    def __getitem__(self, seq_name):
        """Returns the amino acid sequence of a particular sequence name.
        """

        if seq_name not in self.seq_annotations:
            raise AlignmentError('Sequence not found in the alignment: {0}'.format(seq_name))
        else:
            seq_index = self._mapping_r[seq_name]
            seq_as_str = self.sequences[seq_index, ].tostring()
            return Sequence._make((seq_name, seq_as_str, self.seq_annotations[seq_name]))

    def __len__(self):
        """
        Returns the number of sequences in the Alignment.
        """
        return len(self.names)

    # Lookup/Access methods
    @property
    def query(self):
        """Returns the first (query) sequence of the alignment"""
        return self.names[0]

    # Filter Methods
    # def filter(self, ann_type, expression):
    #     """Returns sequences matching a given annotation type & expression

    #     e.g. filter('f_gaps' '>= 0.3')
    #     e.g. filter('OS' '== Homo sapiens')
    #     """

    #     pass

    # Edition methods: annotate/add/remove
    def annotate(self, annotations, overwrite=False):
        """Adds annotations to the alignment

        Args:
            annotations [dict]: Dictionary with per-sequence/generic annotations.
                                If key matches sequence name, added to seq annotation.

            overwrite [bool]: Forces writing of annotations.
                              If False raises an error when annotation key already exists.
        """

        for seq_name in annotations:
            if seq_name in self.seq_annotations:
                _existing_keys = self.seq_annotations[seq_name].viewkeys()
                _new_keys = annotations[seq_name].viewkeys()

                _overlapping_keys = _existing_keys & _new_keys
                if _overlapping_keys and not overwrite:
                    _overlapping_keys = ','.join(_overlapping_keys)
                    emsg = 'Overlapping annotations for sequence {0}: {1}'
                    raise AlignmentError(emsg.format(seq_name, _overlapping_keys))

                self.seq_annotations[seq_name].update(annotations[seq_name])
                self.logger.debug('Annotated seq. {0} with keys {1}'.format(seq_name, _new_keys))
            else:
                raise AlignmentError('Sequence not found in the alignment: {0}'.format(seq_name))

    def add(self, sequence, seq_name, annotations=None):
        """Appends a sequence to the bottom of the alignment

        Args:
            seq_name [str]: Sequence name. Must be unique in the alignment.
            sequence [str]: Aminoacid sequence.
            annotations [dict]: Per-sequence annotations to include.
            overwrite [bool]: forces writing of annotations. If False raises error when
                              annotation key already exists.
        """

        if seq_name not in self.seq_annotations:
            self.seq_annotations[seq_name] = {}
        else:
            raise AlignmentError('Sequence already in the alignment: {0}'.format(seq_name))

        if annotations:
            self.logger.debug('Adding annotations to seq. {0}'.format(seq_name))
            self.annotate(annotations)

        self.sequences = np.vstack((self.sequences, sequence))
        self.logger.debug('New 2D array shape: {0}'.format(self.sequences.shape))

        self.names.append(seq_name)
        self._mapping_r = {seq_i: seq_n for seq_i, seq_n in enumerate(self.names)}

    def remove(self, seq_name):
        """Removes a sequence from the alignment along with its annotations.
        """

        if seq_name not in self.seq_annotations:
            raise AlignmentError('Sequence not found in the alignment: {0}'.format(seq_name))

        seq_index = self._mapping_r[seq_name]
        self.sequences = np.delete(self.sequences, seq_index, axis=0)
        del self.seq_annotations[seq_name]
        self.logger.debug('New 2D array shape: {0}'.format(self.sequences.shape))

        self.names.remove(seq_name)
        self._mapping_r = {seq_i: seq_n for seq_i, seq_n in enumerate(self.names)}

    def trim_query_gaps(self):
        """Removes columns from the alignment that match gaps in the query sequence.
        Modifies alignment in-place."""

        keep_mask = self.sequences[0] != '-'
        self.sequences = self.sequences[:, keep_mask]

        n_gaps = sum(self.sequences[0] == '-')
        self.logger.debug('Removed {0} columns from alignment'.format(n_gaps))

    # Statistics methods
    def query_coverage(self, threshold=0.3, annotate=True):
        """Calculates the coverage of the query.

        Coverage is defined as the fraction of positions (columns) in the alignment with less than
        30% (default) gaps throughout the alignment.

        Also defines an additional high-level annotation with '-' for positions below the coverage
        threshold.

        Args:
            threshold [float]: fraction of gaps maximum to consider the position 'gappy'.
            annotate [bool]: adds a high-level annotation with '-' for positions with low coverage.
        Returns:
            coverage [float]
        """

        gappy_columns = ((self.sequences == '-').sum(0) / self.sequences.shape[0]) >= threshold
        sequence_coverage = 1 - (gappy_columns.sum() / self.sequences.shape[0])

        if annotate:
            self.aln_annotations['coverage'] = ''.join(['-' if p else 'x' for p in gappy_columns])

        return sequence_coverage

    def sequence_gaps(self, annotate=True):
        """Calculates the fraction of gaps in each sequence.

        Args:
            - annotate [bool]: adds the result as an additional per-sequence annotation.
        """

        gaps_per_row = (self.sequences == '-').sum(1) / self.sequences.shape[1]

        if annotate:
            for seq_name, seq_gaps in zip(self.names, gaps_per_row):
                self.seq_annotations[seq_name]['f_gaps'] = seq_gaps

        return gaps_per_row


    def alignment_variability(self, annotate=True):
        """Calculates the per-position variability in the alignment, expressed as the number of
        unique elements in the alignment column divided by the number of possible elements (21).

        Returns a number between 0 (all 21 elements observed) and 1 (only one element observed).

        Args:
            - annotate [bool]: adds the result as an additional alignment annotation.
        """

        def _count_unique(array):
            """Returns the number of unique elements in the array"""
            return len(set(array))

        per_column_variation = 1 - (np.apply_along_axis(_count_unique, 0, self.sequences) / 21)

        if annotate:
            self.aln_annotations['variability'] = per_column_variation

        return per_column_variation.sum() / self.sequences.shape[1]

    # Plotting functions
