#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# dependencies: scikit-bio >= 0.2.3, < 0.3.0

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os import makedirs
from os.path import join, exists
import numpy
import numpy.testing as npt
import pandas as pd

from skbio.util import remove_files
import skbio.io

from distance_method import (preprocess_data,
                             parse_blast,
                             normalize_distances,
                             cluster_distances,
                             detect_outlier_genes,
                             launch_blast,
                             launch_diamond,
                             distance_method)


class DistanceMethodTests(TestCase):
    """ Tests for distance-method HGT detection """

    def setUp(self):
        """ Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        self.target_proteomes_dir = join(self.working_dir, "DB")
        if not exists(self.target_proteomes_dir):
            makedirs(self.target_proteomes_dir)

        # species 1
        self.species_1_fp = join(self.target_proteomes_dir, "species_1.fasta")
        with open(self.species_1_fp, 'w') as tmp:
            tmp.write(species_1)

        # species 2
        self.species_2_fp = join(self.target_proteomes_dir, "species_2.fasta")
        with open(self.species_2_fp, 'w') as tmp:
            tmp.write(species_2)

        # species 3
        self.species_3_fp = join(self.target_proteomes_dir, "species_3.fasta")
        with open(self.species_3_fp, 'w') as tmp:
            tmp.write(species_3)

        # species 4
        self.species_4_fp = join(self.target_proteomes_dir, "species_4.fasta")
        with open(self.species_4_fp, 'w') as tmp:
            tmp.write(species_4)

        # blast alignments (query vs. all species)
        self.blast_fp = join(self.working_dir, "blast.txt")
        with open(self.blast_fp, 'w') as tmp:
            tmp.write(blast_alignments)

        # blast alignments (query vs. all species)
        self.phylip_fp = join(self.working_dir, "distances.txt")
        with open(self.phylip_fp, 'w') as tmp:
            tmp.write(phylip_output)

        # list of files to remove
        self.files_to_remove = [self.species_1_fp,
                                self.species_2_fp,
                                self.species_3_fp,
                                self.species_4_fp,
                                self.blast_fp,
                                self.phylip_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.working_dir)

    def assert_frames_equal(self, actual, expected, use_close=False):
        """
        Custom assert_frames_equal() function for testing pandas DataFrame.

        See
        http://nbviewer.ipython.org/gist/jiffyclub/ac2e7506428d5e1d587b
        for details.

        Compare DataFrame items by index and column and
        raise AssertionError if any item is not equal.

        Ordering is unimportant, items are compared only by label.
        NaN and infinite values are supported.

        Parameters
        ----------
        actual : pandas.DataFrame
        expected : pandas.DataFrame
        use_close : bool, optional
            If True, use numpy.testing.assert_allclose instead of
            numpy.testing.assert_equal.

        """
        if use_close:
            comp = npt.assert_allclose
        else:
            comp = npt.assert_equal

        assert (isinstance(actual, pd.DataFrame) and
                isinstance(expected, pd.DataFrame)), \
            'Inputs must both be pandas DataFrames.'

        for i, exp_row in expected.iterrows():
            assert i in actual.index, 'Expected row {!r} not found.'.format(i)

            act_row = actual.loc[i]

            for j, exp_item in exp_row.iteritems():
                assert j in act_row.index, \
                    'Expected column {!r} not found.'.format(j)

                act_item = act_row[j]

                try:
                    comp(act_item, exp_item)
                except AssertionError as e:
                    raise AssertionError(
                        e.message + '\n\nColumn: {!r}\nRow: {!r}'.format(j, i))

    def test_preprocess_data(self):
        """ Test functionality of preprocess_data()
        """
        gene_map, ref_db, species = preprocess_data(self.working_dir,
                                                    self.target_proteomes_dir,
                                                    ['fa', 'fasta', 'faa'])
        gene_map_exp = {'G1_SE001': '0_0', 'G1_SE002': '1_0',
                        'G1_SE003': '2_0', 'G1_SE004': '3_0',
                        '0_0': 'G1_SE001', '1_0': 'G1_SE002',
                        '2_0': 'G1_SE003', '3_0': 'G1_SE004',
                        'G2_SE001': '0_1', 'G2_SE002': '1_1',
                        'G2_SE003': '2_1', 'G2_SE004': '3_1',
                        '0_1': 'G2_SE001', '1_1': 'G2_SE002',
                        '2_1': 'G2_SE003', '3_1': 'G2_SE004',
                        'G3_SE001': '0_2', 'G3_SE002': '1_2',
                        'G3_SE003': '2_2', 'G3_SE004': '3_2',
                        '0_2': 'G3_SE001', '1_2': 'G3_SE002',
                        '2_2': 'G3_SE003', '3_2': 'G3_SE004',
                        'G4_SE001': '0_3', 'G4_SE002': '1_3',
                        'G4_SE003': '2_3', 'G4_SE004': '3_3',
                        '0_3': 'G4_SE001', '1_3': 'G4_SE002',
                        '2_3': 'G4_SE003', '3_3': 'G4_SE004',
                        'G5_SE001': '0_4', 'G5_SE002': '1_4',
                        'G5_SE003': '2_4', 'G5_SE004': '3_4',
                        '0_4': 'G5_SE001', '1_4': 'G5_SE002',
                        '2_4': 'G5_SE003', '3_4': 'G5_SE004'}
        ref_db_exp = {}
        for seq in skbio.io.read(self.species_1_fp, format='fasta'):
            ref_db_exp[seq.metadata['id']] = seq
        for seq in skbio.io.read(self.species_2_fp, format='fasta'):
            ref_db_exp[seq.metadata['id']] = seq
        for seq in skbio.io.read(self.species_3_fp, format='fasta'):
            ref_db_exp[seq.metadata['id']] = seq
        for seq in skbio.io.read(self.species_4_fp, format='fasta'):
            ref_db_exp[seq.metadata['id']] = seq
        num_species_exp = 4
        self.assertDictEqual(gene_map, gene_map_exp)
        self.assertDictEqual(ref_db, ref_db_exp)
        self.assertEqual(species, num_species_exp)

    def test_parse_blast(self):
        """ Test functionality of parse_blast()
        """
        hits_exp = {'G1_SE001': ['G1_SE001', 'G1_SE002', 'G1_SE003',
                                 'G1_SE004'],
                    'G4_SE001': ['G4_SE001', 'G4_SE002', 'G4_SE003',
                                 'G4_SE004'],
                    'G3_SE001': ['G3_SE001', 'G3_SE002', 'G3_SE003',
                                 'G3_SE004'],
                    'G5_SE001': ['G5_SE001', 'G5_SE002', 'G5_SE003',
                                 'G5_SE004'],
                    'G2_SE001': ['G2_SE001', 'G2_SE002', 'G2_SE003',
                                 'G2_SE004']}
        gene_map = {'G1_SE001': '0_0', 'G1_SE002': '1_0', 'G1_SE003': '2_0',
                    'G1_SE004': '3_0', '0_0': 'G1_SE001', '1_0': 'G1_SE002',
                    '2_0': 'G1_SE003', '3_0': 'G1_SE004', 'G2_SE001': '0_1',
                    'G2_SE002': '1_1', 'G2_SE003': '2_1', 'G2_SE004': '3_1',
                    '0_1': 'G2_SE001', '1_1': 'G2_SE002', '2_1': 'G2_SE003',
                    '3_1': 'G2_SE004', 'G3_SE001': '0_2', 'G3_SE002': '1_2',
                    'G3_SE003': '2_2', 'G3_SE004': '3_2', '0_2': 'G3_SE001',
                    '1_2': 'G3_SE002', '2_2': 'G3_SE003', '3_2': 'G3_SE004',
                    'G4_SE001': '0_3', 'G4_SE002': '1_3', 'G4_SE003': '2_3',
                    'G4_SE004': '3_3', '0_3': 'G4_SE001', '1_3': 'G4_SE002',
                    '2_3': 'G4_SE003', '3_3': 'G4_SE004', 'G5_SE001': '0_4',
                    'G5_SE002': '1_4', 'G5_SE003': '2_4', 'G5_SE004': '3_4',
                    '0_4': 'G5_SE001', '1_4': 'G5_SE002', '2_4': 'G5_SE003',
                    '3_4': 'G5_SE004'}
        hits = {}
        parse_blast(self.blast_fp, hits, gene_map)
        self.assertDictEqual(hits, hits_exp)

    def test_normalize_distances(self):
        """ Test functionality of normalize_distances()

        Phylip alignments (row IDs symbolize species_gene):
        2_1         0.000000  0.379562  0.473355  0.521700
        3_1         0.379562  0.000000  0.587981  0.660393
        0_1         0.473355  0.587981  0.000000  0.722046
        1_1         0.521700  0.660393  0.722046  0.000000

        Z-score normalized by rows:
        [[        nan -1.33276108  0.25673322  1.07602786]
         [-1.36991607         nan  0.38082407  0.989092  ]
         [-1.19162122 -0.06375679         nan  1.25537801]
         [-1.3488877   0.30650842  1.04237928         nan]]

        Re-ordered by rows and columns to correspond to ascending species names
        [[        nan  1.25537801 -1.19162122 -0.06375679]
         [ 1.04237928         nan -1.3488877   0.30650842]
         [ 0.25673322  1.07602786         nan -1.33276108]
         [ 0.38082407  0.989092   -1.36991607         nan]]
    """
        num_species = 4
        i = 0
        species_set_dict = {}
        species_set_dict_exp = {'IIII': 1}
        gene_bitvector_map = {}
        gene_bitvector_map_exp = {0: 'IIII'}
        full_distance_matrix = numpy.zeros(
            shape=(1, num_species, num_species), dtype=float)
        full_distance_matrix_exp = numpy.array(
            [[[numpy.nan, 1.25537801, -1.19162122, -0.06375679],
              [1.04237928, numpy.nan, -1.3488877, 0.30650842],
              [0.25673322, 1.07602786, numpy.nan, -1.33276108],
              [0.38082407, 0.989092, -1.36991607, numpy.nan]]])
        normalize_distances(phylip_fp=self.phylip_fp,
                            full_distance_matrix=full_distance_matrix,
                            num_species=num_species,
                            full_distance_matrix_offset=i,
                            species_set_dict=species_set_dict,
                            gene_bitvector_map=gene_bitvector_map)
        numpy.testing.assert_almost_equal(full_distance_matrix[0][0],
                                          full_distance_matrix_exp[0][0])
        self.assertDictEqual(species_set_dict, species_set_dict_exp)
        self.assertDictEqual(gene_bitvector_map, gene_bitvector_map_exp)

    def test_cluster_distances(self):
        """ Test functionality of cluster_distances()
        """
        species_set_dict = {'IIIIIIII': 100, 'IIOOOIII': 50, 'IIIIIII0': 10,
                            'OIOIIIII': 5, 'IIIOOIII': 8, 'OOOOOIOO': 12}
        gene_clusters_list_exp = [('IIIIIIII', ['IIIIIIII', 'IIIIIII0',
                                                'IIIOOIII', 'OIOIIIII']),
                                  ('IIOOOIII', ['IIOOOIII', 'OOOOOIOO'])]
        gene_clusters_list_act = cluster_distances(
            species_set_dict=species_set_dict, species_set_size=30,
            hamming_distance=2)
        self.assertTrue(len(gene_clusters_list_exp),
                        len(gene_clusters_list_act))
        for core_cluster_exp in gene_clusters_list_exp:
            self.assertTrue(core_cluster_exp in gene_clusters_list_act)
        for core_cluster_act in gene_clusters_list_act:
            self.assertTrue(core_cluster_act in gene_clusters_list_exp)

    def test_detect_outlier_genes(self):
        """ Test functionality of detect_outlier_genes()
        """
        species_set = ['IIII']
        gene_bitvector_map = {0: 'IIII', 1: 'IIII', 2: 'IIII',
                              3: 'IIII', 4: 'IIII'}
        full_distance_matrix = numpy.array(
            [[[numpy.nan, 1.20467207, 0.03920422, -1.24387629],
              [0.70710678, numpy.nan, -1.41421356, 0.70710678],
              [0.70710678, 1.41421356, numpy.nan, -0.70710678],
              [1.24387629, 1.20467207, 0.03920422, numpy.nan]],
             [[numpy.nan, 1.26889551, -1.175214, -0.09368151],
              [1.16820922, numpy.nan, -1.27436935, 0.10616013],
              [0.50122985, 0.89462587, numpy.nan, -1.39585572],
              [0.55177142, 0.85179386, -1.40356529, numpy.nan]],
             [[numpy.nan, 1.33958186, -1.06239803, -0.27718382],
              [1.2373387, numpy.nan, -0.02558867, -1.21175004],
              [0.284687, 1.05732936, numpy.nan, -1.34201637],
              [0.78533243, 0.62588164, -1.41121407, numpy.nan]],
             [[numpy.nan, 1.38826553, -0.92766886, -0.46059667],
              [1.15415521, numpy.nan, -1.2848518, 0.13069659],
              [0.28409367, 1.05773152, numpy.nan, -1.34182519],
              [0.26316662, 1.0717693, -1.33493592, numpy.nan]],
             [[numpy.nan, 1.25537801, -1.19162122, -0.06375679],
              [1.04237928, numpy.nan, -1.3488877, 0.30650842],
              [0.25673322, 1.07602786, numpy.nan, -1.33276108],
              [0.38082407, 0.989092, -1.36991607, numpy.nan]]])
        outlier_genes_exp = set([0])
        outlier_genes = detect_outlier_genes(
            species_set=species_set,
            gene_bitvector_map=gene_bitvector_map,
            full_distance_matrix=full_distance_matrix,
            stdev_offset=1.5,
            outlier_hgt=0.5,
            num_species=4,
            total_genes=5)
        self.assertSetEqual(outlier_genes, outlier_genes_exp)

    def test_launch_blast(self):
        """Test functionality of launch_blast()
        """
        align_exp = [{'qseqid': 'G1_SE001',
                      'sseqid': 'G1_SE002',
                      'pident': 60.92,
                      'length': 888,
                      'mismatch': 328,
                      'gapopen': 5,
                      'qstart': 1,
                      'qend': 870,
                      'sstart': 1,
                      'send': 887,
                      'evalue': 0.000000e+00,
                      'bitscore': 1098,
                      'qcovs': 100},
                     {'qseqid': 'G2_SE001',
                      'sseqid': 'G2_SE002',
                      'pident': 53.64,
                      'length': 494,
                      'mismatch': 229,
                      'gapopen': 0,
                      'qstart': 1,
                      'qend': 494,
                      'sstart': 1,
                      'send': 494,
                      'evalue': 0.000000e+00,
                      'bitscore': 566,
                      'qcovs': 100},
                     {'qseqid': 'G3_SE001',
                      'sseqid': 'G3_SE002',
                      'pident': 64.66,
                      'length': 116,
                      'mismatch': 40,
                      'gapopen': 1,
                      'qstart': 1,
                      'qend': 115,
                      'sstart': 1,
                      'send': 116,
                      'evalue': 2.9999999999999994e-56,
                      'bitscore': 164,
                      'qcovs': 100},
                     {'qseqid': 'G4_SE001',
                      'sseqid': 'G4_SE002',
                      'pident': 48.29,
                      'length': 292,
                      'mismatch': 147,
                      'gapopen': 1,
                      'qstart': 1,
                      'qend': 288,
                      'sstart': 1,
                      'send': 292,
                      'evalue': 1.9999999999999996e-106,
                      'bitscore': 305,
                      'qcovs': 100},
                     {'qseqid': 'G5_SE001',
                      'sseqid': 'G5_SE002',
                      'pident': 50.00,
                      'length': 670,
                      'mismatch': 320,
                      'gapopen': 6,
                      'qstart': 2,
                      'qend': 663,
                      'sstart': 1,
                      'send': 663,
                      'evalue': 0.000000e+00,
                      'bitscore': 674,
                      'qcovs': 99}]
        df_exp = pd.DataFrame(align_exp,
                              columns=['qseqid', 'sseqid', 'pident', 'length',
                                       'mismatch', 'gapopen', 'qstart', 'qend',
                                       'sstart', 'send', 'evalue', 'bitscore',
                                       'qcovs'])
        out_file_fp = launch_blast(self.species_1_fp,
                                   self.species_2_fp,
                                   self.working_dir)
        df_act = skbio.io.read(out_file_fp, format='blast+6',
                               into=pd.DataFrame,
                               columns=['qseqid', 'sseqid', 'pident', 'length',
                                        'mismatch', 'gapopen', 'qstart',
                                        'qend', 'sstart', 'send', 'evalue',
                                        'bitscore', 'qcovs'])
        self.assert_frames_equal(df_exp, df_act)

    def test_launch_diamond(self):
        """Test functionality of launch_diamond()
        """
        align_exp = [{'qseqid': 'G1_SE001',
                      'sseqid': 'G1_SE002',
                      'pident': 60.9,
                      'length': 888,
                      'mismatch': 328,
                      'gapopen': 5,
                      'qstart': 1,
                      'qend': 870,
                      'sstart': 1,
                      'send': 887,
                      'evalue': 0.000000e+00,
                      'bitscore': 1092.8},
                     {'qseqid': 'G2_SE001',
                      'sseqid': 'G2_SE002',
                      'pident': 53.6,
                      'length': 494,
                      'mismatch': 229,
                      'gapopen': 0,
                      'qstart': 1,
                      'qend': 494,
                      'sstart': 1,
                      'send': 494,
                      'evalue': 3.1999999999999995e-160,
                      'bitscore': 550.1},
                     {'qseqid': 'G3_SE001',
                      'sseqid': 'G3_SE002',
                      'pident': 64.7,
                      'length': 116,
                      'mismatch': 40,
                      'gapopen': 1,
                      'qstart': 1,
                      'qend': 115,
                      'sstart': 1,
                      'send': 116,
                      'evalue': 2.2999999999999996e-45,
                      'bitscore': 166.4},
                     {'qseqid': 'G4_SE001',
                      'sseqid': 'G4_SE002',
                      'pident': 48.3,
                      'length': 292,
                      'mismatch': 147,
                      'gapopen': 1,
                      'qstart': 1,
                      'qend': 288,
                      'sstart': 1,
                      'send': 292,
                      'evalue': 4.799999999999999e-84,
                      'bitscore': 296.2},
                     {'qseqid': 'G5_SE001',
                      'sseqid': 'G5_SE002',
                      'pident': 50.1,
                      'length': 669,
                      'mismatch': 319,
                      'gapopen': 6,
                      'qstart': 2,
                      'qend': 662,
                      'sstart': 1,
                      'send': 662,
                      'evalue': 1.4999999999999994e-192,
                      'bitscore': 657.9}]
        df_exp = pd.DataFrame(align_exp,
                              columns=['qseqid', 'sseqid', 'pident', 'length',
                                       'mismatch', 'gapopen', 'qstart', 'qend',
                                       'sstart', 'send', 'evalue', 'bitscore'])
        out_file_fp = launch_diamond(self.species_1_fp,
                                     self.species_2_fp,
                                     self.working_dir,
                                     tmp_dir=self.working_dir)
        df_act = skbio.io.read(out_file_fp, format='blast+6',
                               into=pd.DataFrame,
                               columns=['qseqid', 'sseqid', 'pident',
                                        'length', 'mismatch', 'gapopen',
                                        'qstart', 'qend', 'sstart', 'send',
                                        'evalue', 'bitscore'])
        self.assert_frames_equal(df_exp, df_act)

    def test_distance_method(self):
        """ Test functionality of distance_method_main()
        """
        output_hgt_fp = join(self.working_dir, "hgt_result.txt")
        distance_method(self.species_1_fp,
                        self.target_proteomes_dir,
                        self.working_dir,
                        output_hgt_fp,
                        'diamond')
        hgt_exp = []
        hgt_act = []
        with open(output_hgt_fp, 'r') as output_hgt_f:
            for line in output_hgt_f:
                if line.startswith('#'):
                    continue
                if line not in ['\n', '\r\n']:
                    hgt_act.append(line.strip().split()[0])
        self.assertListEqual(hgt_exp, hgt_act)

    def test_distance_method_pass_alignments(self):
        """ Test functionality of distance_method_main() with alignments
        """
        output_hgt_fp = join(self.working_dir, "hgt_result.txt")
        distance_method(self.species_1_fp,
                        self.target_proteomes_dir,
                        self.working_dir,
                        output_hgt_fp,
                        'diamond',
                        tabular_alignments_fp=self.blast_fp)
        hgt_exp = []
        hgt_act = []
        with open(output_hgt_fp, 'r') as output_hgt_f:
            for line in output_hgt_f:
                if line.startswith('#'):
                    continue
                if line not in ['\n', '\r\n']:
                    hgt_act.append(line.strip().split()[0])
        self.assertListEqual(hgt_exp, hgt_act)


phylip_output = """    4
2_1         0.000000  0.379562  0.473355  0.521700
3_1         0.379562  0.000000  0.587981  0.660393
0_1         0.473355  0.587981  0.000000  0.722046
1_1         0.521700  0.660393  0.722046  0.000000
"""

blast_alignments = """G1_SE001    G1_SE001    100.00  862 0   0   1   862 1  \
 862 0.0  1803   100
G2_SE001   G2_SE001    100.00  494 0   0   1   494 1   494 0.0  1023   100
G3_SE001   G3_SE001    100.00  115 0   0   1   115 1   115 1e-85     239   100
G4_SE001   G4_SE001    100.00  288 0   0   1   288 1   288 0.0   599   100
G5_SE001   G5_SE001    100.00  663 0   0   1   663 1   663 0.0  1377   100
G1_SE001   G1_SE002    58.11   888 345 6   1   862 1   887 0.0  1048   100
G2_SE001   G2_SE002    53.64   494 229 0   1   494 1   494 0.0   566   100
G3_SE001   G3_SE002    64.66   116 40  1   1   115 1   116 3e-56     164   100
G4_SE001   G4_SE002    48.29   292 147 1   1   288 1   292 2e-106    305   100
G5_SE001   G5_SE002    50.00   670 320 6   2   663 1   663 0.0   674   99
G1_SE001   G1_SE003    66.74   869 280 4   1   862 1   867 0.0  1191   100
G2_SE001   G2_SE003    64.65   495 174 1   1   494 1   495 0.0   655   100
G3_SE001   G3_SE003    68.97   116 35  1   1   115 1   116 2e-59     172   100
G4_SE001   G4_SE003    57.40   277 113 2   1   273 1   276 1e-117    335   95
G5_SE001   G5_SE003    58.01   674 262 9   1   663 1   664 0.0   769   100
G1_SE001   G1_SE004    65.25   872 291 6   1   862 1   870 0.0  1142   100
G2_SE001   G2_SE004    59.92   494 198 0   1   494 1   494 0.0   578   100
G3_SE001   G3_SE004    66.12   121 35  2   1   115 1   121 6e-59     171   100
G4_SE001   G4_SE004    55.23   277 120 1   1   273 1   277 7e-111    318   95
G5_SE001   G5_SE004    56.06   685 270 8   1   663 1   676 0.0   738   100
G1_SE001   G1_SE004    100.00  862 0   0   1   862 1   862 0.0  1803   100
G2_SE001   G2_SE004    59.92   494 198 0   1   494 1   494 0.0   578   100
G3_SE001   G3_SE004    66.12   121 35  2   1   115 1   121 6e-59     171   100
G4_SE001   G4_SE004    55.23   277 120 1   1   273 1   277 7e-111    318   95
G5_SE001   G5_SE004    56.06   685 270 8   1   663 1   676 0.0   738   100
"""

species_1 = """>G1_SE001
FDDSSLLEIFTSNNSNSSFSEPTVQLASYAEADPVEAASLSGILGQCTRVRHMMSSVTREVMPLQSTRSAKYVGPGV
PPFATAGQGGGDEQFKMADTPCKGVKMEKLKWAEDRHKPLVFLIGDAMYLMVPAENKITQYYNGICNGAGEVWDHLF
YKAECLHCFGFVGESVAYGNNGWSVADVGTVGTKGAGYMVYESLHATIPYALNGRQTDGLRLTYEPEDGSMLAANAI
PYGCVGPDCGDIGEVQSYGQMSNLGEYHLATFKLERDKMRVSAKDAKDSEYPVDGQEGFTDSSDGKGVDVYGPGQHA
YARLVVGKRDRQHATLAEMAEDGYADKMEPRCAQQPATINYNAGEVVGEERITTDIIAREYMFTKLTWNKTSPGYNY
VGAVQSTLLDFPGLWTATNVSREEQAKIHHPEGNVPDHLFCQPNNPPRDYPAKLILFLGILTSTIKSPAETWDAGLS
GQDSKIEMVKLHPLYHIDSSYAPMLNKHSSCIGCPTPLMLPPSAGKLLMLRPHEGTTTATESESYDTGSSAKFFLCY
SPDPVEIFGVPMMQAHNYHPKSVWFHLGNVLKHLGGSKDTSWRGLIVHMPRLLLEQLDAFTELGNGNHKYDSEISND
LGTEGLVALKRRILAQAYAAPNANDFYIGHDTLCAPFIASRKIWAWGKTQVSLEKGNAWAHAVLSPWIIKKEVAQGT
AITALIKSRPIDLPGNGIIGTHHDRPIGAMMVSAKAEEALAATAALPTTALAVSYETASDARQGLIGGLHSSPQFAP
AITGLINYLIERTDNVDLHMAFYVLHVGIVPQKYLARKSTRRGTCWDMHQGCLNTACRSLPAPNAQYHIPISKTLTL
TTAMHKTCIDLAKVWLGDAGGPL
>G2_SE001
GKLKNSIIDPWPGDPFAIASSDQTAAVAIVHAASEYISHYGYGYKAQLMLKIDIQESCNANGGAGGGGCQYAWWTAW
ASIFTQSPDVSISQTRVIYFTSAALGLIGFWSMLLGGAFRGGAEAWNAKADLKQAKVKTRKSAFFNNKNELASVPDI
VLPYPADKSDSSSDMKYFGSSMSKKMIAGYTPAAARPRVTITVEELKMSAKAEYLEDLLKNSPHLHGEIGDSTEKML
LSMTAQCKCTATSEALSYKDQKGGRDAAAPKKDFHGSCGVTFPYGCCYPEKKTAADEIVNLALGVCSSNLKVLQRPG
NIKVEAYIDACVVLDGNVKTGDGESRITLDEIHPFSVLLGEGNISKKVTGTHIGSHFDSITIPIGGQFGLAGVELIT
YQADSKDAIGRAYDKKPIILWFQGVAHELGGPIVPAADETIRIPDYITFVEFKHFDPSTSVVCEDDAAKLDENDKER
ANETKVQEEHSLKAVPTRKRLGARAKSIPFEL
>G3_SE001
VVEDNNQGAPGVVQIFYGNGTLHQEDCFSGPQAIGPGDASPGTLIQVVRGRKTHTEFVNALIKGTDNAPTGERVHIQ
WGLLMPPNFLGPEVKTNLYLDKNFKCFKQFGSIQVKAS
>G4_SE001
LNMYVATSHQEFTGQLYDGKKPTPLVDSPPMNDCQRMSWLFMHTLNTRYKSNDLANGEVRLKAQKHVYVQSFKAATA
YASKVILIEVVTLEQVKSSTLALANAFEKISVAVYKQLLRYATVSETTPGSVLQVEVGARDGILFDGEMLVHSDEAN
SIWGLVLYKGSAKSKLHFGYLFPVTAVIGKVTFPKFKRHPNAGYVDGGLPALKMAFTLTFKFSSHFYPRVQDQRFKD
WINVFHVPYFWGDVKKQRALNLGSTLELLNGVVSDPCEYRLLEETGLGGKAKNAVRT
>G5_SE001
DTGFASEDYVGVEWHYTEVIVVLESPKDRRYKSTAFPKEIGCGYGTLAENRSIWERGRTEPNANELISSSPPLVFPP
IMAHAHGNLPPYGQSKWGSVTWYKLLLALIAETYAVLNVLGLTADPLPRRLGVGVTNVCGAHVLFHDEEKSKEGQTT
VKSSLLIDIEGAALKSILQYFLTEASDNGKDTKENLKPCVYEIRWMSEEGSIAVPDYLSATDFACGSGFVFLMMVNK
FGYFEIRDGQICGGEGTLLILVLGQEIPDEKYFMAKGTTFGRNDFESNAMNHQKVVMTILPNNWPISGVKTDTADQI
LDGCFGFIPLPSARKATFAVDSALGTGSHLIKRTGTNAMVIYVVIVLVCATLVPGSPNYLTGIMLSDVQLLVCDALS
DSVKWLFAPTKLLEIKPTYIVTADMHTSTETKAQKDVIGKRMDFASNGLNTRAIKIEFQLLYSMAYGFLGFLVLRAC
TCGKFAQIVNDVCATLVWGDALGYLSATNQNTLEYGTTGHSENTESFELNPYKIENQQTDEAPRIANKIALKRNGAA
AGRMAANYTLGDFYLLTEYSCNNCKVTDGAVFDYATGVERGLDRHTEVQLVTPEPLLTGEAQNKHQLAVRVGWLAYA
QFMAPPIEADVTQTSLLPASVTRGYERSGETGSGFTKTALGNEAGAQ
"""

species_2 = """>G1_SE002
FDDVSLLNIVTSDNSQTKFNVPTVQLSGVLEAEPTELAVLSAILAMCNAVRQILPGVTRDPTDLAKRRTARYVGPGV
PPFATAGQGGGDEKVDMNETPCKGFPIGQLIWAEERDKALVFLNSENVYLLVPSENKTEYLEGICNGATNVWGHLFY
RSDCLHCFALVGDSVPSGNGGWQVTDLGAVGTRATGYMIYEHLQAGIPYALDGLQTAGLRITYAPQVANMLPANAIP
YDCIGPDCGEIGEVLAYGHCSSLGEYHLESFKLERDKLKVSAKSATECEYPVGGTIQFTDTSSASGVDVYGPGHHVY
ERLVVAEKDHQLASFAELADDGHADKVDEACAEGSATIEYSKGEEKGEDVIPTTTLFAKTYMTKTVQRGKTSPGYNY
AGTVKSTLLAMPGMWVADNIAYEEQAKIHHPQGNVPNHLFCNPNNLPRDYPAKLILFLSILTAEIKSPMAVWTAGLS
GQDNRIRLIKLEALWHIDVHCAAPKPLPAPSMYLPCLKKSSIPLGVPDELMLPKATGKLLMLRTSHEKSKIATEFIA
YDSNMSVKFFLCYAPDQVEIYGEPFVQSQKSTPKAVWFHLGAFLKHLAGSKDKAWKGLIVHMPRILIEKLECFYDLG
NGNHNYDSEISNDLGTEGLVTLARRILAEAYAAPNDSEFYIGQDTLCSPFVASRKAWSWGRIQVSLEKGDEHKNAVL
SPWRIKKEVAMGAPITSLKKSGPIGLPANGIVGSHHDKPIKERVVSANATEALAALGALARTQIAASTQIASQEREG
KYGALQMTPKFAYSITGSISYPIERTDKVDILMARYVLHEGIIEQSYLARGSANRGLCWVYGQGCLNTACSPTGPTL
PAPNAQYHLPLSSTLVTESAMGKTCIDPVKGWPGDAGGPL
>G2_SE002
GKLKDSFFHPWPTDPLATSRTSQTSALAIIEMPSSYITHYGYGFWFHKMLKYNFSDSCQAKGGAGDGGCNRDWWTSW
ADCFQQNPDISVAQSRVMYVESAALGLMGFWFLHLGGAFKGGEEKWHVKTDLRCSNVAPKSAGFLKDKQRLASFPDI
RTPKGKDKSDTSTNMNYFGSLLSQRMISGYTPKAGRPRITLIVEELKMAKKVKYLTDFLKMIPHLHHANGDHEEKEL
LGLTPQAKCSATNQFSCFQHSKTGTDPAALKMAFNGSCGFTLSYQACFEYHASAGSAILHITVQICPSGLKVLKRQG
AVKVAASAEYAVVLDGLSKYYDGSSRIIADEHKPLQVLLGAGTLNNGVMGTHWGSHVESMTIAIGGEFGLASVELSS
YQADSCNALGSAYDNKPIILWFQAVAHNLGGPMIPSNNATVRIPQYVSFVEFEHFNPSTGVVCQDDQNRLDQDDRER
GDREGVQESHGLKALPTLKRLSKQALSIAFDL
>G3_SE002
VVEDNSQGAPGVVSVFFGNGTLHQADCFLGPQAVGGGDTHPGTLVQVIRGRKAHTELVNALVGTSNNAPTGDQVHIQ
WGISVTSEPMGPEASNKLFKIDKNFRCFKQFGSAEVWTS
>G4_SE002
MNMYIATADMQFTGKFHDGEKSAPALDAPALKEDEMFQWIFTHKKSTKYDSADLLKGEVQLKDRKHVYIDDFTADAA
FSSKVFVLEVASSKQTKTAALGLRMVLEKSNIAVVKERLRYATVYDKTRGTLLIVGVGYRNGVLFDGEMLVHNEEGN
AIWGLVLAKGDASTLMHVGYIVASASVIVSVTFRKFNRKPNDGYSDGGLPTLKASFSWTFRFCNHFWTEIFDQRLTR
QVQDIIAVISEPFYYSDAKRERKLHLGNTLKVLSGVVSDPCEFKLMDESKLAGAQKTLART
>G5_SE002
TGVESDSYIGEDPHYSNIVTVIDNPNEQKYKTNAFPQSMGCGYGSVAESHKIFGRGRHEPHLKHVMYSSPPFIFPPI
LSTANGNEPAYGQSQRGSAAWYKWLSKVSNVKIPLISSSSAALLVLGLTSLAVPKGLGSGWTAGCGNFVIFHGEEGN
EEGQATVASDLLIGVERGVLRAVLGYYLTESSDNKKDTEEDLRPCIYIIRWPSREGLQSVRNIVKATESAVTPGIVF
FMMISKFATFILGDGKVCGGAGMLLINILGGEIPEVKYFQAKGTTFGYGGFATGSMDHDNTVTVPPNNWPMTGGLTS
IAFQLLSDTFGFVPKPIAQRATFALEEELGTNAQSIKTTDTKALVIYVVGLIVCASLTPGQSHVLHGIILSDVRLVV
CDAASRGVQYAPTRELEIKPTYIFPTDSHESMATQAKTELLGAQGEFAANGLGIEKEHQSDYSFAYGFLGFISFRAG
TCGKFAEIVTDTCKSMVFGAQLRTLWASKETTLEYATDGHTAQSESWPLGPFKVEKRSTDEAKSVATNVGLKVNGED
AARETEEYRLGDFYLLTEYVDNCKVTEGKWLDYASSVEKGSDRHMKVQMIAPKPLISGKGRAGSQASNRVGWLEYRN
HMASPLESEVTRSHLDGACVRRGYDRVGEMGSGLTKNSLNISEAAAH
"""

species_3 = """>G1_SE003
FDDISLLEIVTSDNSNSSLSVATVQLMTYAEALPVEAASLSGQLAQCTTVRHIVSSVSRDVMPLQSTRSARYVGPEV
PPFKTAGQGGGDEEFNMSESPCKGLPMDKLKWAEERHKALVFLLGDAMFLLVPAENKTKWYKGICNGSGEVWDHLFY
KSECLHCFALVGMSVAYGNNGLQVAHVGTVGTKGAGYMIFEWLGAAIPYALNGLQTDGLRITYEPQIGAMLLANIPY
GCIGPDCGDVGEVHNYGQCSNLREYHLATKKLERDKMRVTAKDAKECEYPVEGQEGFTDSSDGSGVDVFGPGSHAFA
RLVVGERDHQHATLAELAEDGYADKMENTCAQGPGTINYNAGEEVGEEVIPTDIIARTYMLKKLQRNKTSPGYNYVG
TVRSTLLAMPGLWTATEVASDEQAKIHHPQGNVPDHLLCKPNSPKQDYPAKLIVFLGILTSTITSPAKVWDAGLSGQ
EDRIELIKLEPLYHINIDYAPLLQRHSICVGCPTPLMVEKSAGKLLMLRTPHEGSTTATEFESYDTGSSVKFFLCTS
PDAIEIYGVPLLQSHDSHPKSVWFHLGNFLKHLGGTKDTSWKGLIVHMPRFVLEKLEAFTQLGKGNHKYDSEISNDL
GTEGLIALTRLILAEASAAPQLNDVYIGQDTLCKPFIASRKAWAWGKIQISLEKGNAWKHSVLSPWIIKKEVAQGAA
VTAKTKSGPIDLPGNGIVGPHKDKPIDCRMVSAHAEEALAALAALPRTTYAVSTETAREAREGLYGALHMSPQFAPA
ITGLLNYIIERTDDVDLLMAFYVLHVGIIEQLYLARKSTNRGLCWDFGQGCLNHACLPAPNAQFHVPISKTLTLETA
LHKTCIDLLKGWLGDSSGAL
>G2_SE003
GKLKDSFFDPWPGDPLATARSHQTAALDIVHAASSYITHYGYGYKSHKMLKYDVSDSCEAKGGAGGGGCQHDGWTTW
AAIFTQSPDISIAQNRVIYFETAALGMMGFWFMELGGVFRGGEEAWTVRAELKCANGAPKKATFFTDKQKLASVPDI
VTPKGVSKFDSSTNMKYMGSVMSKKGLSGYIPAAGRPRITMEVEELKWASKAEYLEDFLKMSPVLHGEIGDHTKKEL
IGLTPQCKCEKTSQALSFKDQKSGVDPAALKKAFGHQSCGISFSYGACFEEKKAGADQIVHLSVQVCSNYLKVLKKQ
HAIKVEAYFESCVELDGLSKYYDGRSPVIIDEHRPLSVLLGAGQLTRKVTGTHIGTHISSMTIAIGGQFGLASVELW
NYQADACDAIGSAYDSEPIILWFQGAAHPLGGPMIPANDETVRVADYVTFVEFEHYDPHTGIVTEDDAAKLDQSERD
RGDETKVQEAHSLEAVPTLKRLGAQAKSVPFEL
>G3_SE003
VVEDNNQGAPGVVYVFHGNGTLHDEDCFTSPRAIGPPDSHPGTLVKVIHGRKAHTEFVNALIGGSENAPTGDRVHIQ
RGIVMPPEPLGPEGQANLDGLDKNFRCLKQLGTAEVRAS
>G4_SE003
LNMYVAKADVQFTGQLHDGKKPSPLFSAPALNDDERLEWVFTHIQSTRYNHGDLQKGQVQLKVRKHVYIQEFEADSA
FSAKIIVLEVKSLESKSAALALKAAFELSNVAVYKQRLRYAKVYAETRGVILLVQVGAREGILFDGAMLIHDDEADV
IVGLALHKGSAQSIMHIGLMFPAAAVIGKVTYCKFTRKPNDGYVNGGLPALKMAFPLSFKFSDHFFDEIQDQAFIRF
VKDWIAIVQVSYFYGDIQKQRKTHLGSTIELLSGVVSDPCEFKFLYFVLLHLRLAITFSEKQLGSLDLLSKMPTLID
ETNLAGADKNLAKT
>G5_SE003
NTGAESVPYVGEEWKYTNVVAILENPESQYYKSSAYPKEIGCAYGGRAPSHSIWERGRHEPNADHLLASSSPLIFAP
VAAHARGDEPPYGQSKWGTLTIHKWLAKTDSKLKLLLLSATYETFCVLGLTAAAVPKDLGKGVTNVCGNVVLFHDAK
GSVEGECTIKSSLLIGIEGAALREILEYFLKQASDKKKDTNTNLKPCLYIIRWASDEGLISVTNYLKATESAIVCGF
VFFMMISKLIYFVFKDGKVCGGEGMVLILILNGEIEEVKYFVAKGTTFGRNGFATNCLSYQGVMTAIPNEWPITGGL
TNTAFQLLDGSFGGVPLPTARKATFANDTSRGTNAQIIKTTDTRSELVIYVVKLIVCSTMNPGQPNLLHGIMLSDNR
LVVCNAASSGVKYTPTKLLEIKPSYYLPGDSHQSTKTKSQKEVIGLRLEFAANGLSIEVEFALSYSFGYGFLGILSL
RAANCGKFAAIVNDTCKTNTWRDQLNVIWAPAEGTLQLATTGQSENTESFDLGPYKVENKHSDEAPRIASKIGLKTN
GIDAGREPEEYAIRGDYYLLTDYCSNLEVTYGKVFDFAAGTEKGLDRHMEVQLITPSPMLSGKGKDQSVRAGWLAYA
NFMAPPLKREVTQSSLDGASVTRGYECHGEMGLSLTTTSLGITERGAQ
"""

species_4 = """>G1_SE004
FDDSSLLEIFTSNNSNSSFSEPTVQLASYAEADPVEAASLSGILGQCTRVRHMMSSVTREVMPLQSTRSAKYVGPGV
PPFATAGQGGGDEQFKMADTPCKGVKMEKLKWAEDRHKPLVFLIGDAMYLMVPAENKITQYYNGICNGAGEVWDHLF
YKAECLHCFGFVGESVAYGNNGWSVADVGTVGTKGAGYMVYESLHATIPYALNGRQTDGLRLTYEPEDGSMLAANAI
PYGCVGPDCGDIGEVQSYGQMSNLGEYHLATFKLERDKMRVSAKDAKDSEYPVDGQEGFTDSSDGKGVDVYGPGQHA
YARLVVGKRDRQHATLAEMAEDGYADKMEPRCAQQPATINYNAGEVVGEERITTDIIAREYMFTKLTWNKTSPGYNY
VGAVQSTLLDFPGLWTATNVSREEQAKIHHPEGNVPDHLFCQPNNPPRDYPAKLILFLGILTSTIKSPAETWDAGLS
GQDSKIEMVKLHPLYHIDSSYAPMLNKHSSCIGCPTPLMLPPSAGKLLMLRPHEGTTTATESESYDTGSSAKFFLCY
SPDPVEIFGVPMMQAHNYHPKSVWFHLGNVLKHLGGSKDTSWRGLIVHMPRLLLEQLDAFTELGNGNHKYDSEISND
LGTEGLVALKRRILAQAYAAPNANDFYIGHDTLCAPFIASRKIWAWGKTQVSLEKGNAWAHAVLSPWIIKKEVAQGT
AITALIKSRPIDLPGNGIIGTHHDRPIGAMMVSAKAEEALAATAALPTTALAVSYETASDARQGLIGGLHSSPQFAP
AITGLINYLIERTDNVDLHMAFYVLHVGIVPQKYLARKSTRRGTCWDMHQGCLNTACRSLPAPNAQYHIPISKTLTL
TTAMHKTCIDLAKVWLGDAGGPL
>G2_SE004
GKLKDSFFDPRPGNPLAVARSHQKAAHAIIHAASNSITHYGYGFRSHKMLKYDVSDAVEAKGGAGGGGCQHDRWKTW
AEIFTQSSDISIAQSRVVYFESAALGLMGFWFMHLGGAFHGGEKAWNVKADLKCANVAPKKATIGRWEVPREGVEDY
FTDEKKQAVFSSTDMKYMGSVMAKKGMSGYTPKAGSRRICVSVEEIKMASDAEYLEEFLKQSPHLRGEIGDHTKKQL
VGMNPQCNCSRTSKALSFKEQSGGVDPAELKKPFHGSCGVTFSYGACFEEKKSGADEIVELQVQTCVSHLKVLKREG
ALKVEAYIESCVVLDGLSKYYDGRSCVVVDEHRPLRVLLGEGTKSKKVTGTHIGTHIKSMTIAIGGQFGLASVESTS
YTADVCDAIGSAYDSKPIVLWFEGAEKEQGGTVIPSNDETARLPDYVTFVEFKHFDPSTGLLCDGVAAKLDQDEKER
GSETKVQEGHALSAVSTLKHLGAQAKSIPFEL
>G3_SE004
VVEDNNQGAPGVVSVFYGDSQLHQEDCFTSPKAIGPGDRHPGTLVQVIRGRKVHTEFTNALIGATEGARTGDRVHIQ
WGIIMEGNAAPPEPLGPESNSNLFGVDKNFRCFKQLGSAEARAS
>G4_SE004
MNMYVATSDIQFSGQVHEGKRPSPLFDSPALGDDKRLRCVVTHIQSTRWDSGDLQKGEVHLKARKHSYIQTLEADSS
FSAKVTVLDIESLSQSKSAALALRARLEKGDVAVYAQRLKSATVYAKTRGTLLLVEVGARMGILYDGEFLIHADEAN
AAVGLVLHKGSAQSVMPIGYLFPPAAVIGKVTFCKFTRAPNDGYVDGELPALKMSFALSFKFTSYFFPEVQDQSFNR
FVKNWIAIVQVAYFYGDIQRQRRTHLPDTIELLSGVVIDPCEGHLLYFVLLHLRLDITYVESRMGSLKQLTKMSTLI
DESNLSGSEKNLTST
>G5_SE004
DTGVESVEYVGEEWHYTTVVPDLENPEKENYKSSTFPKSIGCGFGNLAQSHTIGERGRHEPNQDHLLAQSPPLIFPP
VLAHAKGNEPPYGQSQWGVATWYKWLAKATSKLKLPLIASTYVLLLVLRFSALALPKDLGKGLANVCGNAVLFHDAK
GREENQAIIITSLLITIEGQALREIKDYFLTSAPTNQKDTDTNLQPCVYVMRRASEEGLISVTNYLKATASAIVCGF
VFFMMLAAFVYFKCESGKVCGREGMLQILILGREIWDGKYFLAKKTTFGTNGFATQCMDYQRVLTPLPNEQAAAPGF
CGKDHSWPINAGRTNTAFQLLDGVFGFRSLPKARKATFAVDTALGTNAELIKTADTKSLVIYLVKLITCATLMPGQP
NLLHGIMLCDTRLVVCDAASSGVKFAPCKLLEIKPNYMLPADSHTSTKTKSQKEVLGLRVEFADNGLSIQVEFALSY
SFGYGFLGFLSLRAGHCGTFAEAVNDACETIVWRDRKGVVWATNENSCAVTGHSTDTESFKLGTYNVENNHTDEAPR
IANKIGLKVTGIAAGREPEEYSIGDFYLVTNGCNNVQVIHGKIFDFASGIESGLDRHMEAQLVTPKCLLTGAGKAQV
QLSIRVGWLMYANLMAPPLKEDVHQDSLDGVSVTRGYECSGEMGLGLTETSMGITDAGAH
"""


if __name__ == '__main__':
    main()
