# ----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join

from skbio.util import remove_files

from benchmark.compute_accuracy import (parse_expected_transfers,
                                        parse_observed_transfers,
                                        compute_accuracy)


class ComputeAccuracyTests(TestCase):
    """ Tests for compute_accuracy.py """

    def setUp(self):
        """ Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        # ALF logfile
        self.alf_log_fp = join(self.working_dir, "alf_log.txt")
        with open(self.alf_log_fp, 'w') as tmp:
            tmp.write(alf_log_file)

        # HGT summary
        self.hgt_sum_fp = join(self.working_dir, "hgt_summary.txt")
        with open(self.hgt_sum_fp, 'w') as tmp:
            tmp.write(hgt_summary)

        # list of files to remove
        self.files_to_remove = [self.alf_log_fp,
                                self.hgt_sum_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.working_dir)

    def test_parse_expected_transfers(self):
        """ Test functionality of parse_expected_transfers()
        """
        exp_trans = [('2', '166', '1', '4122'), ('3', '3015', '2', '4123'),
                     ('3', '1318', '5', '4124'), ('5', '2736', '1', '4125'),
                     ('2', '2594', '6', '4126'), ('1', '349', '6', '4127'),
                     ('2', '2286', '1', '4128'), ('3', '1481', '1', '4129'),
                     ('10', '4027', '2', '4160'), ('9', '4001', '10', '4161'),
                     ('3', '2409', '1', '4162')]
        with open(self.alf_log_fp, 'U') as alf_f:
            observed_transfers = parse_expected_transfers(
                ground_truth_f=alf_f)
        self.assertItemsEqual(observed_transfers, exp_trans)

    def test_parse_observed_transfers(self):
        """ Test functionality of parse_observed_transfers()
        """
        exp_hgts = {'Jane 4': ['1000', '1105'],
                    'Consel': ['1000', '1105'],
                    'RIATA-HGT': ['1000', '1105'],
                    'RANGER-DTL': ['1000', '1105'],
                    'T-REX': ['1000', '1105']}
        with open(self.hgt_sum_fp, 'U') as hgt_f:
            obs_hgts = parse_observed_transfers(observed_hgts_f=hgt_f,
                                                pvalue_cutoff=0.05)
        self.assertDictEqual(obs_hgts, exp_hgts)

    def test_compute_accuracy(self):
        """ Test functionality of compute_accuracy()
        """
        ground_truth_trans =\
            [('2', '166', '1', '4122'), ('3', '1481', '1', '4129'),
             ('10', '4027', '2', '4160')]
        tools_trans = {'Jane 4': ['166', '4027'],
                       'Consel': ['166', '4027'],
                       'RIATA-HGT': ['166', '4027'],
                       'RANGER-DTL': ['166', '4027'],
                       'T-REX': ['166', '4027']}
        tools_accuracy_exp =\
            {'Jane 4': (2, 0, 1, 1.0, 0.6666666666666666, 0.8),
             'Consel': (2, 0, 1, 1.0, 0.6666666666666666, 0.8),
             'RIATA-HGT': (2, 0, 1, 1.0, 0.6666666666666666, 0.8),
             'RANGER-DTL': (2, 0, 1, 1.0, 0.6666666666666666, 0.8),
             'T-REX': (2, 0, 1, 1.0, 0.6666666666666666, 0.8)}
        tools_accuracy_obs = compute_accuracy(ground_truth_trans, tools_trans)
        self.assertDictEqual(tools_accuracy_obs, tools_accuracy_exp)


alf_log_file = """Synthetic evolution
-------------------

first organism (biological sequences): 4120 proteins with 314 aa (average)
height of tree: 9.28, Speciations: 10
LGT rate: 0.00


time 1.1111: deletion of length 1 in org/gene 1/1929
something interesting is happening here...
time 1.5940: speciation event of organism 1 to organism 2
time 1.6708: lgt from organism 2 with gene 166 to organism 1, now gene 4122
    orthologues replacement, gene 166 in organism 1 deleted
time 1.8385: speciation event of organism 1 to organism 3
time 1.9703: insertion of length 1 in org/gene 1/439
time 2.5678: speciation event of organism 1 to organism 4
time 2.7886: lgt from organism 3 with gene 3015 to organism 2, now gene 4123
    orthologues replacement, gene 3015 in organism 2 deleted
time 2.8812: insertion of length 2 in org/gene 1/2172
time 2.9393: speciation event of organism 1 to organism 5
time 2.9576: lgt from organism 3 with gene 1318 to organism 5, now gene 4124
    orthologues replacement, gene 1318 in organism 5 deleted
time 2.9825: lgt from organism 5 with gene 2736 to organism 1, now gene 4125
    orthologues replacement, gene 2736 in organism 1 deleted
time 2.9783: deletion of length 1 in org/gene 1/2860
time 3.1447: speciation event of organism 1 to organism 6
time 3.2272: speciation event of organism 2 to organism 7
time 3.2440: lgt from organism 2 with gene 2594 to organism 6, now gene 4126
    orthologues replacement, gene 2594 in organism 6 deleted
time 3.3550: speciation event of organism 1 to organism 8
time 3.4516: lgt from organism 1 with gene 349 to organism 6, now gene 4127
    orthologues replacement, gene 349 in organism 6 deleted
time 3.5742: lgt from organism 2 with gene 2286 to organism 1, now gene 4128
    orthologues replacement, gene 2286 in organism 1 deleted
time 3.7624: lgt from organism 3 with gene 1481 to organism 1, now gene 4129
    orthologues replacement, gene 1481 in organism 1 deleted
time 6.7650: deletion of length 1 in org/gene 1/43
time 4.9727: insertion of length 1 in org/gene 1/105
time 3.4662: deletion of length 2 in org/gene 1/480
time 4.5764: deletion of length 2 in org/gene 1/697
time 6.2927: deletion of length 1 in org/gene 1/906
time 5.9451: insertion of length 4 in org/gene 1/1874
time 5.8455: deletion of length 1 in org/gene 1/3283
time 7.1311: speciation event of organism 1 to organism 10
time 7.2242: lgt from organism 10 with gene 4027 to organism 2, now gene 4160
    orthologues replacement, gene 4027 in organism 2 deleted
time 7.3389: lgt from organism 9 with gene 4001 to organism 10, now gene 4161
    orthologues replacement, gene 4001 in organism 10 deleted
time 7.4675: lgt from organism 3 with gene 2409 to organism 1, now gene 4162
    orthologues replacement, gene 2409 in organism 1 deleted
time 8.7236: deletion of length 2 in org/gene 1/1831
time 7.4542: insertion of length 2 in org/gene 1/3433
time 3.5660: deletion of length 2 in org/gene 2/2592
time 4.3849: deletion of length 1 in org/gene 3/3957
time 5.8385: insertion of length 4 in org/gene 4/579


41200 genes in 10 species generated
62 lgt events resulting in 96 genes
"""

hgt_summary = """#number of HGTs detected
#\tgene ID\tT-REX\tRANGER-DTL\tRIATA-HGT\tJane 4\tConsel
0\t1000\t1\t1\t1\t1\t0.99 0.01
1\t1001\t0\t0\t0\t0\t0.00 0.00
2\t1105\t1\t1\t1\t1\t0.99 0.01
"""


if __name__ == '__main__':
    main()
