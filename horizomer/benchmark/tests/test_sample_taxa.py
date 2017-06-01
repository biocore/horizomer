#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join, dirname, realpath

from horizomer.benchmark.sample_taxa import sample_taxa


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

    def tearDown(self):
        # there isn't any file to remove at the moment
        # but in the future there will be
        rmtree(self.working_dir)

    def test_sample_taxa(self):
        """ Test sample taxa from hit table
        """
        obs = sample_taxa(self.hit_table_fp,
                          self.prot2tax_dict_fp)
        exp = set()
        with open(self.exp_taxa_fp, 'r') as f:
            for line in f:
                exp.add(line.rstrip())
        self.assertSetEqual(obs, exp)


if __name__ == '__main__':
    main()
