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
from os import makedirs, listdir
from os.path import join, dirname, realpath
from skbio import io
from click.testing import CliRunner

from horizomer.benchmark.sample_genes import (
    sample_genes,
    write_genes,
    _main)


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
        self.ref_faa_dir = join(datadir, 'ref_faa')
        self.write_genes_list = join(datadir, 'CrDC.write_genes.list')

        # intermediate files
        self.out_fa_dir = join(self.working_dir, 'out_fa')
        makedirs(self.out_fa_dir)
        self.out_genes_fp = join(self.working_dir, 'genes.list')

    def tearDown(self):
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

    def test_write_genes(self):
        genes = sample_genes(self.ortho_groups_fp, min_taxa_cutoff=4.0)
        write_genes(genes, self.ref_faa_dir, self.out_fa_dir,
                    self.out_genes_fp)
        # test number of output FASTA files
        obs = sorted(listdir(self.out_fa_dir))
        exp = sorted(['%s.fa' % ogid for ogid in genes])
        self.assertListEqual(obs, exp)
        # test output gene list content
        with open(self.out_genes_fp, 'r') as f:
            obs = f.read()
        with open(self.write_genes_list, 'r') as f:
            exp = f.read()
        self.assertEqual(obs, exp)
        # test FASTA file content
        for seq in io.read(join(self.out_fa_dir, 'OG0000017.fa'),
                           format='fasta'):
            exp = 'GCF_000160655.1|WP_040356123.1'
            self.assertEqual(seq.metadata['id'], exp)
            exp = ('MYRKHYAADVTETLDGQTVQVAGWVHRRRDHGGVIFIDLRDRSGLVQIVIDPDTADAF'
                   'ALAEQVRNEYCLAIEGRVRLRPAGTENPDLASGKIEILGKQLTVLSKSEPLPFQLDED'
                   'NVSEEIRLKHRTIDLRRDVMQKNLILRSKVAASLRRYLDEHGFMDIETPMLTKATPEG'
                   'ARDYLVPSRTHPGKFFALPQSPQLFKQMLMMSGFDRYYQIVRCFRDEDLRADRQPEFT'
                   'QLDIETSFLEEEDILQIMEPMIRGIFKEHLGVELANPFPRMTYREAMRRYASDKPDLR'
                   'IPLELVDIDDLVKNSGFKVFASVAAQDNGRVVALKIPGGAKLTRKEIDDYTAYVARYG'
                   'AKGLAYIKVNDATNVEGLQSPIVKFLTTEGGAEGAIALDIIKRVDAQNGDLIFFGADK'
                   'ASIVNDAIGALRIKVGHDLNMLTCDWAPLWVVDFPMFEYDEKDGRWYSMHHPFTQPKT'
                   'ANLDELDTNPGDVLSRAYDMVLNGTEIGGGSIRIHRDDMQQRVFKSLGIGAEEAQEKF'
                   'GFLLNALKYGCPPHGGIAFGLDRLIMLMAGAKSIRDVMAFPKTQTAWCPLTDAPSEAS'
                   'EAQLRELHIRKRQVEKSE')
            self.assertEqual(str(seq), exp)
            break

    def test__main(self):
        params = ['--ortho-groups-fp', self.ortho_groups_fp,
                  '--min-taxa-cutoff', 4.0,
                  '--input-faa-dir', self.ref_faa_dir,
                  '--output-fa-dir', self.out_fa_dir,
                  '--output-genes-fp', self.out_genes_fp]
        res = CliRunner().invoke(_main, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output, 'Number of gene families sampled: 6.\n')


if __name__ == '__main__':
    main()
