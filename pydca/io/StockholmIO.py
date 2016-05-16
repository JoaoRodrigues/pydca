# Copyright 2016 by Joao Rodrigues, Jason Wang, Jessica Zhao.  All rights reserved.
# This code is part of the pydca distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Classes to read/write Stockholm alignment files.
"""

import logging
import os
import re

from ..core.Alignment import Alignment

class ParserError(Exception):
    """Custom Parser Exception"""
    pass

class StockholmReader(object):
    """
    Class to parse multiple sequence alignment files in Stockholm format into an Alignment
    object.
    """

    def __init__(self, hmmer_fmt=True):

        self.logger = logging.getLogger(__name__)

        self._seq_data = {} # ID: sequence & GS/GR records
        self._aln_data = {'query_id': None} # GF record data
        self._ordered_seq_ids = [] # list of IDs as they appear in the sequence section
        self.is_hmmer = hmmer_fmt

    def read(self, stockholm_file):
        """Public method to parse file. Returns Alignment object"""

        fpath = os.path.abspath(stockholm_file)
        if not os.path.isfile(fpath):
            raise IOError("File not found: {0}".format(fpath))

        msa = self._parse(fpath)
        return msa

    # Private line-parsing methods

    # Per-sequence annotations
    def _parse_GS_record(self, line):
        """Private method to parse GS records"""

        fields = line.split()
        seqname = fields[1]
        feat_name = fields[2].upper()
        feat_data = " ".join(fields[3:])

        if seqname not in self._seq_data:
            self._seq_data[seqname] = {}

        if feat_name in self._seq_data[seqname]:
            emsg = 'Duplicated feature {0} for sequence {1}'.format(feat_name, seqname)
            raise ParserError(emsg)
        else:
            self._seq_data[seqname][feat_name] = feat_data

        self.logger.debug('Parsed GS: {0} => {1}'.format(feat_name, feat_data))

    def _parse_GS_record_hmmer(self, line):
        """Private method to parse HMMER-specific GS records"""

        fields = line.split()
        seqname = fields[1]
        feat_name = fields[2].upper()
        feat_data = " ".join(fields[3:])

        if seqname not in self._seq_data:
            self._seq_data[seqname] = {}

        # Parse DE line for additional annotations
        # Split at XX= occurences
        if feat_name == 'DE':
            nested_annotations = re.split('([A-Z][A-Z])=', feat_data)
            nested_annotations.insert(0, 'DE') # re-insert 'DE' identifier

            # iterate over annotations
            num_annotations = len(nested_annotations)

            assert num_annotations % 2 == 0, "Unexpected format at line: {0}".format(line)

            for index in range(0, num_annotations, 2):
                feat_name = nested_annotations[index]
                feat_data = nested_annotations[index + 1]
                if feat_name in self._seq_data[seqname]:
                    emsg = 'Duplicated feature {0} for sequence {1}'.format(feat_name, seqname)
                    raise ParserError(emsg)
                else:
                    self.logger.debug('Parsed nested GS: {0} = {1}'.format(feat_name, feat_data))
                    self._seq_data[seqname][feat_name] = feat_data

        else:
            if feat_name in self._seq_data[seqname]:
                emsg = 'Duplicated feature {0} for sequence {1}'.format(feat_name, seqname)
                raise ParserError(emsg)
            else:
                self._seq_data[seqname][feat_name] = feat_data

            self.logger.debug('Parsed GS: {0} => {1}'.format(feat_name, feat_data))

    def _parse_sequence(self, line):
        """Private method to parse sequence data"""

        fields = line.split()
        seqname = fields[0]
        seqdata = fields[1]

        # Record first sequence as query in metadata
        if not self._aln_data['query_id']:
            self._aln_data['query_id'] = seqname

        if seqname not in self._seq_data:
            self._seq_data[seqname] = {}

        if 'seq' not in self._seq_data[seqname]:
            self._seq_data[seqname]['seq'] = list(seqdata)
            self._ordered_seq_ids.append(seqname)
            self.logger.debug('Parsed Sequence: {0}'.format(seqname))
        else:
            self._seq_data[seqname]['seq'].extend(seqdata)
            self.logger.debug('Extended Sequence: {0}'.format(seqname))

    # Per-Column Annotations
    def _parse_GR_record(self, line):
        """Private method to parse GR annotations"""

        fields = line.split()
        seqname = fields[1]
        feat_name = fields[2].upper()
        feat_data = " ".join(fields[3:])

        if seqname not in self._seq_data:
            raise ParserError('Orphan GR/GC record: {0}'.format(line))

        elif feat_name in self._seq_data[seqname]:
            # Duplicates are possible if the sequence data is wrapped.
            self._seq_data[seqname][feat_name] += feat_data

        else:
            self._seq_data[seqname][feat_name] = feat_data

    # Alignment / Consensus Annotations
    def _parse_GF_record(self, line):
        """Private method to parse GF records"""

        fields = line.split()
        feat_name = fields[1].upper()
        feat_data = " ".join(fields[2:])

        if feat_name not in self._aln_data:
            self._aln_data[feat_name] = feat_data
        else:
            # Append to existing feature
            self._aln_data[feat_name] += feat_data

        self.logger.debug('Parsed GF: {0} => {1}'.format(feat_name, feat_data))

    def _parse_GC_record(self, line):
        """Private method to parse GC annotations"""

        fields = line.split()
        feat_name = fields[1].upper()
        feat_data = " ".join(fields[2:])

        if feat_name in self._aln_data:
            # Duplicates are possible if the sequence data is wrapped.
            self._aln_data[feat_name] += feat_data
        else:
            self._aln_data[feat_name] = feat_data

    def _parse(self, fpath):
        """
        High-level parser function.

        Raises:
            ParserError: an error occured while parsing the Stockholm file
        """

        # Define GS parser based on user input
        if self.is_hmmer:
            _GS_line_parser = self._parse_GS_record_hmmer
        else:
            _GS_line_parser = self._parse_GS_record

        # Actually parse the file
        with open(fpath, 'rU') as handle:

            # Verify presence of Stockholm header on first line
            line = handle.next()
            if not line.startswith('# STOCKHOLM 1.0'):
                raise ParserError('Missing Stockholm file header')

            for line in handle:

                line = line.strip()

                # Parse features & Sequences
                if line[0:4] == '#=GF':
                    self._parse_GF_record(line)

                elif line[0:4] == '#=GS':
                    _GS_line_parser(line)

                elif line[0:4] == '#=GR':
                    self._parse_GR_record(line)

                elif line[0:4] == '#=GC':
                    self._parse_GC_record(line)

                elif line[0:2] == '//':
                    break

                elif line:
                    self._parse_sequence(line)

        # Build Alignment object
        # Use order from sequence data, not GS records (i.e. self._ordered_seq_ids)
        # Ensure consistency, raise warnings for:
        #   - sequences without annotation
        #   - annotations for missing sequences (skip)
        _sequences = []
        for seqname in self._ordered_seq_ids:
            seq = self._seq_data[seqname]['seq']
            _sequences.append(seq)
            del self._seq_data[seqname]['seq']

        n_sequences = len(_sequences)
        self.logger.debug('Creating Alignment with {0} sequences'.format(n_sequences))

        msa = Alignment(_sequences, self._ordered_seq_ids,
                        seq_annotations=self._seq_data,
                        aln_annotations=self._aln_data)

        return msa

# class StockholmWriter():

#     """
#     Wrapper class that converts multiple sequence sequences from matrix form to Stockholm
#     """

#     def __init__(self, MSA, stockholm_filename):

#         # Setup logging with module name
#         self.gr_flags = ["SS", "SA", "TM", "PP", "LI", "AS", "pAS", "sAS", "IN"]
#         self.logger = logging.getLogger(__name__)
#         self.output_alignment(MSA, stockholm_filename)

#     def _output_GF(self, MSA, output_file):

#         for key, value in MSA.metadata.items():
#             output_file.write("#=GF {0} {1}\n".format(key, value))
#         output_file.write("\n")

#     def _output_GS(self, MSA, output_file):

#         for annotation in MSA.annotations:
#             # "ID" key must be in the annotations dictionary
#             assert("ID" in annotation)

#             for key, value in annotation.items():
#                 if key not in (["ID"] + self.gr_flags):
#                     output_file.write("# GS {0} {1} {2}\n".format(annotation['ID'], key, value))
#         output_file.write("\n")

#     def _output_sequence_and_GR(self, MSA, output_file):

#         line = 200 # limit number of characters per line
#         num_sequence_sections = len(MSA.sequences)/line #integer division
#         if num_sequence_sections == 0:
#             num_sequence_sections += 1

#         seq_struct = "{0:<35s} {1}\n"
#         gr_struct = "#=GR {0:<30s} {1}\n"

#         for sect in range(int(num_sequence_sections)):
#             for i, sequence in enumerate(MSA.sequences):
#                 sequence = "".join(sequence)

#                 end_index = min((sect+1)*line, MSA.sequences.shape[1]-sect*line)
#                 seq = seq_struct.format(MSA.annotations[i]['ID'], sequence[sect*line:end_index])
#                 output_file.write(seq)

#                 for gr_entry in self.gr_flags:
#                     if gr_entry in MSA.annotations[i]:
#                         gr_header = "{0} {1}".format(MSA.annotations[i]['ID'], gr_entry)
#                         output_file.write((gr_struct.format(gr_header, MSA.annotations[i][gr_entry][sect*line:end_index])))

#             # Separate sequence sections with a space
#             if sect != num_sequence_sections-1: # No space following last section
#                 output_file.write("\n")

#     def output_alignment(self, MSA, stockholm_filename):
#         """
#         Outputs the MSA object into a Stockholm-formatted file

#         Args:
#             MSA: MSA object to be outputted
#             stockholm_filename: an absolute or relative path to desired output file name

#         Raises:
#             ParseError: an error occured while parsing the Stockholm file
#         """

#         output_file = open(stockholm_filename, 'w')

#         # Write Header
#         output_file.write("# STOCKHOLM 1.0\n")

#         # Write per file flags
#         self._output_GF(MSA, output_file)

#         # Write per sequence flags
#         self._output_GS(MSA, output_file)

#         # Write sequences and per residue flags
#         self._output_sequence_and_GR(MSA, output_file)

#         # End Stockholm
#         output_file.write("//\n")
