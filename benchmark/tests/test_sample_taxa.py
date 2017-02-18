# -----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os import remove
from os.path import join, dirname, realpath

from benchmark.sample_taxa import sample_taxa


class SampleTaxaTests(TestCase):
    """ Test for sample_taxa.py """

    def setUp(self):
        """ Set up working directory and test files
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()

        # test data directory
        datadir = join(dirname(realpath(__file__)), 'data', 'simplified')

        # test data files
        self.hit_table_fp = join(datadir, 'CrDC.m8')
        self.prot2tax_dict_fp = join(datadir, 'prot2gcf.txt')
        self.exp_taxa_fp = join(datadir, 'CrDC.m8.taxa')
        self.output_taxa_fp = join(self.working_dir, 'output.taxa')
        self.files_to_remove = [self.output_taxa_fp]

    def tearDown(self):
        for file in self.files_to_remove:
            remove(file)
        rmtree(self.working_dir)

    def test_sample_taxa(self):
        """ Test sample taxa from hit table
        """
        sample_taxa(self.hit_table_fp,
                    self.prot2tax_dict_fp,
                    self.output_taxa_fp)
        with open(self.output_taxa_fp, 'r') as f:
            obs = f.read()
        with open(self.exp_taxa_fp, 'r') as f:
            exp = f.read()
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
