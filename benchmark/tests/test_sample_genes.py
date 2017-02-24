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
from os.path import join, dirname, realpath

from benchmark.sample_genes import sample_genes


class SampleGenesTests(TestCase):
    """ Test for sample_taxa.py """

    def setUp(self):
        """ Set up working directory and test files
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()

        # test data directory
        datadir = join(dirname(realpath(__file__)), 'data', 'simplified')

        # test data files
        self.ortho_groups_fp = join(datadir, 'CrDC.orthogroups.tsv')
        self.sample_genes_return = join(datadir, 'CrDC.sample_genes.return')

    def tearDown(self):
        # there isn't any file to remove at the moment
        # but in the future there will be
        rmtree(self.working_dir)

    def test_sample_genes(self):
        """ Test sample taxa from hit table
        """
        with open(self.sample_genes_return, 'r') as f:
            exp = eval(f.read())
        # cutoff as a number
        obs = sample_genes(self.ortho_groups_fp, min_taxa_cutoff=4.0)
        self.assertDictEqual(obs, exp)
        # cutoff as a fraction (5 * 0.8 = 4)
        obs = sample_genes(self.ortho_groups_fp, min_taxa_cutoff=0.8)
        # invalid cutoff
        self.assertRaises(ValueError, sample_genes, self.ortho_groups_fp,
                          min_taxa_cutoff=-1)


if __name__ == '__main__':
    main()
