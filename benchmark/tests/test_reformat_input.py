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
from tempfile import mkstemp, mkdtemp
from os import close, remove
from os.path import join

from skbio import TreeNode, TabularMSA, Protein

from benchmark.reformat_input import (join_trees,
                                      trim_gene_tree_leaves,
                                      species_gene_mapping,
                                      remove_branch_lengths,
                                      id_mapper,
                                      reformat_rangerdtl,
                                      reformat_trex,
                                      reformat_riatahgt,
                                      reformat_jane4,
                                      reformat_treepuzzle,
                                      reformat_genemark,
                                      reformat_egid)


class ReformatInputTests(TestCase):
    """ Test WGS-HGT input reformatting functions """

    def setUp(self):
        """
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()

        # species tree
        f, self.species_tree_fp = mkstemp(prefix='species_',
                                          suffix='.nwk')
        close(f)
        with open(self.species_tree_fp, 'w') as t:
            t.write(species_tree)

        # species tree 2
        f, self.species_tree_2_fp = mkstemp(prefix='species_2_',
                                            suffix='.nwk')
        close(f)
        with open(self.species_tree_2_fp, 'w') as t:
            t.write(species_tree_2)

        # gene tree 1
        f, self.gene_tree_1_fp = mkstemp(prefix='gene_tree_1_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_1_fp, 'w') as t:
            t.write(gene_tree_1)

        # gene tree 2
        f, self.gene_tree_2_fp = mkstemp(prefix='gene_tree_2_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_2_fp, 'w') as t:
            t.write(gene_tree_2)

        # gene tree 3
        f, self.gene_tree_3_fp = mkstemp(prefix='gene_tree_3_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_3_fp, 'w') as t:
            t.write(gene_tree_3)

        # MSA FASTA 3
        f, self.msa_fa_3_fp = mkstemp(prefix='msa_3_',
                                      suffix='.fa')
        close(f)
        with open(self.msa_fa_3_fp, 'w') as t:
            t.write(msa_fa_3)

        # species genome
        f, self.species_genome_fp = mkstemp(prefix='species_',
                                            suffix='.gbk')
        close(f)
        with open(self.species_genome_fp, 'w') as t:
            t.write(species_genome)

        self.files_to_remove = [self.species_tree_fp,
                                self.species_tree_2_fp,
                                self.gene_tree_1_fp,
                                self.gene_tree_2_fp,
                                self.gene_tree_3_fp,
                                self.msa_fa_3_fp,
                                self.species_genome_fp]

    def tearDown(self):
        for file in self.files_to_remove:
            remove(file)
        rmtree(self.working_dir)

    def test_join_trees(self):
        """ Test concatenate Newick trees into one file (species, gene)
        """
        self.output_file = join(self.working_dir, 'output_file.nwk')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        join_trees(gene_tree_1, species_tree, self.output_file)
        with open(self.output_file, 'r') as out_f:
            species_gene_tree_1_obs = out_f.read()
        self.assertEqual(species_gene_tree_1_obs, species_gene_tree_1_exp)

    def test_trim_gene_tree_leaves(self):
        """ Test remove '_GENENAME' from tree leaf names (if exists)
        """
        leaves_exp = ["SE001", "SE002", "SE003", "SE004", "SE005", "SE006",
                      "SE007", "SE008", "SE009", "SE010"]
        gene_tree_2 = TreeNode.read(self.gene_tree_2_fp, format='newick')
        trim_gene_tree_leaves(gene_tree_2)
        leaves_obs = sorted([node.name for node in gene_tree_2.tips()])
        self.assertListEqual(leaves_obs, leaves_exp)

    def test_species_gene_mapping(self):
        """ Test finding the association between species and gene tree leaves
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        mapping_exp = {"SE001": ["SE001_01623", "SE001_04123"],
                       "SE002": ["SE002_01623"], "SE003": ["SE003_01623"],
                       "SE004": ["SE004_01623"],
                       "SE005": ["SE005_01623", "SE005_04123"],
                       "SE006": ["SE006_01623", "SE006_04123"],
                       "SE007": ["SE007_01623"],
                       "SE008": ["SE008_01623", "SE008_04123"],
                       "SE009": ["SE009_01623", "SE009_04123"],
                       "SE010": ["SE010_01623", "SE010_04123"]}
        mapping_obs = species_gene_mapping(gene_tree_1, species_tree)
        self.assertDictEqual(dict(mapping_obs), mapping_exp)

    def test_species_gene_mapping_check_species_labels(self):
        """ Test ValueError raised
        """
        species_tree = TreeNode.read(self.species_tree_2_fp, format='newick')
        gene_tree_3 = TreeNode.read(self.gene_tree_3_fp, format='newick')
        self.assertRaises(ValueError,
                          species_gene_mapping,
                          gene_tree=gene_tree_3,
                          species_tree=species_tree)

    def test_remove_branch_lengths(self):
        """ Test removing branch lengths from tree
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        remove_branch_lengths(tree=species_tree)
        species_tree_exp = (
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),"
            "SE004),SE003),(SE002,SE007));")
        self.assertEqual(species_tree_exp, str(species_tree)[:-1])

    def test_id_mapper(self):
        """ Test id_mapper() functionality
        """
        index = [u'SE001/02297', u'SE002/02297', u'SE003/02297',
                 u'SE004/02297', u'SE005/02297', u'SE006/02297',
                 u'SE008/02297', u'SE009/02297', u'SE010/02297']
        mapping = id_mapper(index)
        mapping_exp = {u'SE008/02297': u'SE008',
                       u'SE003/02297': u'SE003',
                       u'SE009/02297': u'SE009',
                       u'SE005/02297': u'SE005',
                       u'SE001/02297': u'SE001',
                       u'SE006/02297': u'SE006',
                       u'SE004/02297': u'SE004',
                       u'SE010/02297': u'SE010',
                       u'SE002/02297': u'SE002'}
        self.assertDictEqual(mapping_exp, mapping)

    def test_reformat_rangerdtl(self):
        """ Test functionality of reformat_rangerdtl()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        reformat_rangerdtl(gene_tree_1,
                           species_tree,
                           output_tree_fp)
        joined_tree_exp = [
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),"
            "SE004),SE003),(SE002,SE007));\n",
            "(((((((SE001_01623,SE010_01623),SE008_01623),(SE006_01623,"
            "SE009_01623)),SE005_01623),SE004_01623),SE003_01623),"
            "((SE002_01623,SE007_01623),((((SE001_04123,SE010_04123),"
            "SE008_04123),(SE006_04123,SE009_04123)),SE005_04123)));\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            joined_tree_act = output_tree_f.readlines()
        self.assertListEqual(joined_tree_exp, joined_tree_act)

    def test_reformat_trex(self):
        """ Test functionality of reformat_trex()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        reformat_trex(gene_tree_1,
                      species_tree,
                      output_tree_fp)
        reformat_tree_exp = [
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374):1.594016;\n",
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,((SE002:"
            "6.0534057,SE007:6.0534057):0.4589905,((((SE001:"
            "2.1494876,SE010:2.1494876):3.7761166,SE008:"
            "5.9256042):0.2102448,(SE006:5.2329068,SE009:"
            "5.2329068):0.9029422):0.2054233,SE005:6.3412723):"
            "0.1711239):1.174147):1.594016;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_riatahgt(self):
        """ Test functionality of reformat_riatahgt()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nex")
        reformat_riatahgt(gene_tree_1,
                          species_tree,
                          output_tree_fp)
        reformat_tree_exp = [
            "#NEXUS\n", "BEGIN TREES;\n",
            "Tree speciesTree = "
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374):1.594016;\n",
            "Tree geneTree = "
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,((SE002:"
            "6.0534057,SE007:6.0534057):0.4589905,((((SE001:"
            "2.1494876,SE010:2.1494876):3.7761166,SE008:"
            "5.9256042):0.2102448,(SE006:5.2329068,SE009:"
            "5.2329068):0.9029422):0.2054233,SE005:6.3412723):"
            "0.1711239):1.174147):1.594016;\n",
            "END;\n",
            "BEGIN PHYLONET;\n",
            "RIATAHGT speciesTree {geneTree};\n",
            "END;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_jane4(self):
        """ Test functionality of reformat_jane4()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nex")
        reformat_jane4(gene_tree_1,
                       species_tree,
                       output_tree_fp)
        reformat_tree_exp = [
            "#NEXUS\n", "begin host;\n",
            "tree host = "
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),SE004),SE003),"
            "(SE002,SE007));\n", "\n",
            "endblock;\n", "begin parasite;\n",
            "tree parasite = "
            "(((((((SE001_01623,SE010_01623),SE008_01623),(SE006_01623,"
            "SE009_01623)),SE005_01623),SE004_01623),SE003_01623),"
            "((SE002_01623,SE007_01623),((((SE001_04123,SE010_04123),"
            "SE008_04123),(SE006_04123,SE009_04123)),SE005_04123)));\n", "\n",
            "endblock;\n",
            "begin distribution;\n",
            "Range SE010_01623:SE010, SE010_04123:SE010, SE009_01623:SE009, "
            "SE009_04123:SE009, SE008_01623:SE008, SE008_04123:SE008, "
            "SE007_01623:SE007, SE006_01623:SE006, SE006_04123:SE006, "
            "SE005_01623:SE005, SE005_04123:SE005, SE004_01623:SE004, "
            "SE003_01623:SE003, SE002_01623:SE002, SE001_01623:SE001, "
            "SE001_04123:SE001;\n",
            "endblock;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_treepuzzle(self):
        """ Test functionality of reformat_treepuzzle()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_3 = TreeNode.read(self.gene_tree_3_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        output_msa_phy_fp = join(self.working_dir, "gene_tree_3.phy")
        reformat_treepuzzle(gene_tree_3,
                            species_tree,
                            self.msa_fa_3_fp,
                            output_tree_fp,
                            output_msa_phy_fp)
        reformat_tree_exp = [
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374);\n",
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,SE002:"
            "7.6865432);\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)
        msa_fa = TabularMSA.read(output_msa_phy_fp, constructor=Protein)
        labels_exp = [u'SE001', u'SE002', u'SE003', u'SE004', u'SE005',
                      u'SE006', u'SE008', u'SE009', u'SE010']
        labels_act = list(msa_fa.index)
        self.assertListEqual(labels_exp, labels_act)

    def test_reformat_genemark(self):
        """ Test functionality of reformat_genemark()
        """
        reformat_genemark(self.species_genome_fp,
                          self.working_dir)
        reformat_gm_exp = {'fna': [
            ">locus001\n",
            "GATCCTCCATATACAACGGTATCTCCACCTCAGGTTTAGATCTCAACAACGGAACCATTGCCGA"
            "CATGAGACAGTTAGGTATCGTCGAGAGTTACAAGCTAAAACGAGCAGTAGTCAGCTCTGCATCT"
            "GAAGCCGCTGAAGTTCTACTAAGGGTGGATAACATCATCCGTGCAAGACCAAGAACCGCCAATA"
            "GACAACATATGTAACATATTTAGGATATACCTCGAAAATAATAAACCGCCACACTGTCATTATT"
            "ATAATTAGAAACAGAACGCAAAAATTATCCACTATATAATTCAAAGACGCGAAAAAAAAAGAAC"
            "AACGCGTCATAGAACTTTTGGCAATTCGCGTCACAAATAAATTTTGGCAACTTATGTTTCCTCT"
            "TCGAGCAGTACTCGAGCCCTGTCTCAAGAATGTAATAATACCCATCGTAGGTATGGTTAAAGAT"
            "AGCATCTCCACAACCTCAAAGCTCCTTGCCGAGAGTCGCCCTCCTTTGTCGAGTAATTTTCACT"
            "TTTCATATGAGAACTTATTTTCTTATTCTTTACTCTCACATCCTGTAGTGATTGACACTGCAAC"
            "AGCCACCATCACTAGAAGAACAGAACAATTACTTAATAGAAAAATTATATCTTCCTCGAAACGA"
            "TTTCCTGCTTCCAACATCTACGTATATCAAGAAGCATTCACTTACCATGACACAGCTTCAGATT"
            "TCATTATTGCTGACAGCTACTATATCACTACTCCATCTAGTAGTGGCCACGCCCTATGAGGCAT"
            "ATCCTATCGGAAAACAATACCCCCCAGTGGCAAGAGTCAATGAATCGTTTACATTTCAAATTTC"
            "CAATGATACCTATAAATCGTCTGTAGACAAGACAGCTCAAATAACATACAATTGCTTCGACTTA"
            "CCGAGCTGGCTTTCGTTTGACTCTAGTTCTAGAACGTTCTCAGGTGAACCTTCTTCTGACTTAC"
            "TATCTGATGCGAACACCACGTTGTATTTCAATGTAATACTCGAGGGTACGGACTCTGCCGACAG"
            "CACGTCTTTGAACAATACATACCAATTTGTTGTTACAAACCGTCCATCCATCTCGCTATCGTCA"
            "GATTTCAATCTATTGGCGTTGTTAAAAAACTATGGTTATACTAACGGCAAAAACGCTCTGAAAC"
            "TAGATCCTAATGAAGTCTTCAACGTGACTTTTGACCGTTCAATGTTCACTAACGAAGAATCCAT"
            "TGTGTCGTATTACGGACGTTCTCAGTTGTATAATGCGCCGTTACCCAATTGGCTGTTCTTCGAT"
            "TCTGGCGAGTTGAAGTTTACTGGGACGGCACCGGTGATAAACTCGGCGATTGCTCCAGAAACAA"
            "GCTACAGTTTTGTCATCATCGCTACAGACATTGAAGGATTTTCTGCCGTTGAGGTAGAATTCGA"
            "ATTAGTCATCGGGGCTCACCAGTTAACTACCTCTATTCAAAATAGTTTGATAATCAACGTTACT"
            "GACACAGGTAACGTTTCATATGACTTACCTCTAAACTATGTTTATCTCGATGACGATCCTATTT"
            "CTTCTGATAAATTGGGTTCTATAAACTTATTGGATGCTCCAGACTGGGTGGCATTAGATAATGC"
            "TACCATTTCCGGGTCTGTCCCAGATGAATTACTCGGTAAGAACTCCAATCCTGCCAATTTTTCT"
            "GTGTCCATTTATGATACTTATGGTGATGTGATTTATTTCAACTTCGAAGTTGTCTCCACAACGG"
            "ATTTGTTTGCCATTAGTTCTCTTCCCAATATTAACGCTACAAGGGGTGAATGGTTCTCCTACTA"
            "TTTTTTGCCTTCTCAGTTTACAGACTACGTGAATACAAACGTTTCATTAGAGTTTACTAATTCA"
            "AGCCAAGACCATGACTGGGTGAAATTCCAATCATCTAATTTAACATTAGCTGGAGAAGTGCCCA"
            "AGAATTTCGACAAGCTTTCATTAGGTTTGAAAGCGAACCAAGGTTCACAATCTCAAGAGCTATA"
            "TTTTAACATCATTGGCATGGATTCAAAGATAACTCACTCAAACCACAGTGCGAATGCAACGTCC"
            "ACAAGAAGTTCTCACCACTCCACCTCAACAAGTTCTTACACATCTTCTACTTACACTGCAAAAA"
            "TTTCTTCTACCTCCGCTGCTGCTACTTCTTCTGCTCCAGCAGCGCTGCCAGCAGCCAATAAAAC"
            "TTCATCTCACAATAAAAAAGCAGTAGCAATTGCGTGCGGTGTTGCTATCCCATTAGGCGTTATC"
            "CTAGTAGCTCTCATTTGCTTCCTAATATTCTGGAGACGCAGAAGGGAAAATCCAGACGATGAAA"
            "ACTTACCGCATGCTATTAGTGGACCTGATTTGAATAATCCTGCAAATAAACCAAATCAAGAAAA"
            "CGCTACACCTTTGAACAACCCCTTTGATGATGATGCTTCCTCGTACGATGATACTTCAATAGCA"
            "AGAAGATTGGCTGCTTTGAACACTTTGAAATTGGATAACCACTCTGCCACTGAATCTGATATTT"
            "CCAGCGTGGATGAAAAGAGAGATTCTCTATCAGGTATGAATACATACAATGATCAGTTCCAATC"
            "CCAAAGTAAAGAAGAATTATTAGCAAAACCCCCAGTACAGCCTCCAGAGAGCCCGTTCTTTGAC"
            "CCACAGAATAGGTCTTCTTCTGTGTATATGGATAGTGAACCAGCAGTAAATAAATCCTGGCGAT"
            "ATACTGGCAACCTGTCACCAGTCTCTGATATTGTCAGAGACAGTTACGGATCACAAAAAACTGT"
            "TGATACAGAAAAACTTTTCGATTTAGAAGCACCAGAGAAGGAAAAACGTACGTCAAGGGATGTC"
            "ACTATGTCTTCACTGGACCCTTGGAACAGCAATATTAGCCCTTCTCCCGTAAGAAAATCAGTAA"
            "CACCATCACCATATAACGTAACGAAGCATCGTAACCGCCACTTACAAAATATTCAAGACTCTCA"
            "AAGCGGTAAAAACGGAATCACTCCCACAACAATGTCAACTTCATCTTCTGACGATTTTGTTCCG"
            "GTTAAAGATGGTGAAAATTTTTGCTGGGTCCATAGCATGGAACCAGACAGAAGACCAAGTAAGA"
            "AAAGGTTAGTAGATTTTTCAAATAAGAGTAATGTCAATGTTGGTCAAGTTAAGGACATTCACGG"
            "ACGCATCCCAGAAATGCTGTGATTATACGCAACGATATTTTGCTTAATTTTATTTTCCTGTTTT"
            "ATTTTTTATTAGTGGTTTACAGATACCCTATATTTTATTTAGTTTTTATACTTAGAGACATTTA"
            "ATTTTAATTCCATTCTTCAAATTTCATTTTTGCACTTAAAACAAAGATCCAAAAATGCTCTCGC"
            "CCTCTTCATATTGAGAATACACTCCATTCAAAATTTTGTCGTCACCGCTGATTAATTTTTCACT"
            "AAACTGATGAATAATCAAAGGCCCCACGTCAGAACCGACTAAAGAAGTGAGTTTTATTTTAGGA"
            "GGTTGAAAACCATTATTGTCTGGTAAATTTTCATCTTCTTGACATTTAACCCAGTTTGAATCCC"
            "TTTCAATTTCTGCTTTTTCCTCCAAACTATCGACCCTCCTGTTTCTGTCCAACTTATGTCCTAG"
            "TTCCAATTCGATCGCATTAATAACTGCTTCAAATGTTATTGTGTCATCGTTGACTTTAGGTAAT"
            "TTCTCCAAATGCATAATCAAACTATTTAAGGAAGATCGGAATTCGTCGAACACTTCAGTTTCCG"
            "TAATGATCTGATCGTCTTTATCCACATGTTGTAATTCACTAAAATCTAAAACGTATTTTTCAAT"
            "GCATAAATCGTTCTTTTTATTAATAATGCAGATGGAAAATCTGTAAACGTGCGTTAATTTAGAA"
            "AGAACATCCAGTATAAGTTCTTCTATATAGTCAATTAAAGCAGGATGCCTATTAATGGGAACGA"
            "ACTGCGGCAAGTTGAATGACTGGTAAGTAGTGTAGTCGAATGACTGAGGTGGGTATACATTTCT"
            "ATAAAATAAAATCAAATTAATGTAGCATTTTAAGTATACCCTCAGCCACTTCTCTACCCATCTA"
            "TTCATAAAGCTGACGCAACGATTACTATTTTTTTTTTCTTCTTGGATCTCAGTCGTCGCAAAAA"
            "CGTATACCTTCTTTTTCCGACCTTTTTTTTAGCTTTCTGGAAAAGTTTATATTAGTTAAACAGG"
            "GTCTAGTCTTAGTGTGAAAGCTAGTGGTTTCGATTGACTGATATTAAGAAAGTGGAAATTAAAT"
            "TAGTAGTGTAGACGTATATGCATATGTATTTCTCGCCTGTTTATGTTTCTACGTACTTTTGATT"
            "TATAGCAAGGGGAAAAGAAATACATACTATTTTTTGGTAAAGGTGAAAGCATAATGTAAAAGCT"
            "AGAATAAAATGGACGAAATAAAGAGAGGCTTAGTTCATCTTTTTTCCAAAAAGCACCCAATGAT"
            "AATAACTAAAATGAAAAGGATTTGCCATCTGTCAGCAACATCAGTTGTGTGAGCAATAATAAAA"
            "TCATCACCTCCGTTGCCTTTAGCGCGTTTGTCGTTTGTATCTTCCGTAATTTTAGTCTTATCAA"
            "TGGGAATCATAAATTTTCCAATGAATTAGCAATTTCGTCCAATTCTTTTTGAGCTTCTTCATAT"
            "TTGCTTTGGAATTCTTCGCACTTCTTTTCCCATTCATCTCTTTCTTCTTCCAAAGCAACGATCC"
            "TTCTACCCATTTGCTCAGAGTTCAAATCGGCCTCTTTCAGTTTATCCATTGCTTCCTTCAGTTT"
            "GGCTTCACTGTCTTCTAGCTGTTGTTCTAGATCCTGGTTTTTCTTGGTGTAGTTCTCATTATTA"
            "GATCTCAAGTTATTGGAGTCTTCAGCCAATTGCTTTGTATCAGACAATTGACTCTCTAACTTCT"
            "CCACTTCACTGTCGAGTTGCTCGTTTTTAGCGGACAAAGATTTAATCTCGTTTTCTTTTTCAGT"
            "GTTAGATTGCTCTAATTCTTTGAGCTGTTCTCTCAGCTCCTCATATTTTTCTTGCCATGACTCA"
            "GATTCTAATTTTAAGCTATTCAATTTCTCTTTGATC\n"
        ], 'gbk': [
            "LOCUS       locus001   5028 bp   DNA   circular   CON   01-JAN-1"
            "900\n",
            "FEATURES             Location/Qualifiers\n",
            "     gene            1..206\n",
            "                     /locus_tag=gene1\n",
            "     CDS             1..206\n",
            "                     /locus_tag=gene1\n",
            "                     /protein_id=AAA98665.1\n",
            "                     /translation=SSIYNGISTSGLDLNNGTIADMRQLGIVES"
            "YKLKRAVVSSASEAAEVLLRVDNIIRARPRTANRQHM\n",
            "     gene            687..3158\n",
            "                     /locus_tag=gene2\n",
            "     CDS             687..3158\n",
            "                     /locus_tag=gene2\n",
            "                     /protein_id=AAA98666.1\n",
            "                     /translation=MTQLQISLLLTATISLLHLVVATPYEAYPI"
            "GKQYPPVARVNESFTFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSD"
            "ANTTLYFNVILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDP"
            "NEVFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPETSYS"
            "FVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYVYLDDDPISSD"
            "KLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYGDVIYFNFEVVSTTDLF"
            "AISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQDHDWVKFQSSNLTLAGEVPKNF"
            "DKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSANATSTRSSHHSTSTSSYTSSTYTAKISS"
            "TSAAATSSAPAALPAANKTSSHNKKAVAIACGVAIPLGVILVALICFLIFWRRRRENPDDENLP"
            "HAISGPDLNNPANKPNQENATPLNNPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSV"
            "DEKRDSLSGMNTYNDQFQSQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTG"
            "NLSPVSDIVRDSYGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPS"
            "PYNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRL"
            "VDFSNKSNVNVGQVKDIHGRIPEML\n",
            "     gene            complement(3300..4037)\n",
            "                     /locus_tag=gene3\n",
            "     CDS             complement(3300..4037)\n",
            "                     /locus_tag=gene3\n",
            "                     /protein_id=AAA98667.1\n",
            "                     /translation=MNRWVEKWLRVYLKCYINLILFYRNVYPPQ"
            "SFDYTTYQSFNLPQFVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLD"
            "FSELQHVDKDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDR"
            "NRRVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEKLISG"
            "DDKILNGVYSQYEEGESIFGSLF\n",
            "ORIGIN\n",
            "        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac"
            " ggaaccattg\n",
            "       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca"
            " gtagtcagct\n",
            "      121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt"
            " gcaagaccaa\n",
            "      181 gaaccgccaa tagacaacat atgtaacata tttaggatat acctcgaaaa"
            " taataaaccg\n",
            "      241 ccacactgtc attattataa ttagaaacag aacgcaaaaa ttatccacta"
            " tataattcaa\n",
            "      301 agacgcgaaa aaaaaagaac aacgcgtcat agaacttttg gcaattcgcg"
            " tcacaaataa\n",
            "      361 attttggcaa cttatgtttc ctcttcgagc agtactcgag ccctgtctca"
            " agaatgtaat\n",
            "      421 aatacccatc gtaggtatgg ttaaagatag catctccaca acctcaaagc"
            " tccttgccga\n",
            "      481 gagtcgccct cctttgtcga gtaattttca cttttcatat gagaacttat"
            " tttcttattc\n",
            "      541 tttactctca catcctgtag tgattgacac tgcaacagcc accatcacta"
            " gaagaacaga\n",
            "      601 acaattactt aatagaaaaa ttatatcttc ctcgaaacga tttcctgctt"
            " ccaacatcta\n",
            "      661 cgtatatcaa gaagcattca cttaccatga cacagcttca gatttcatta"
            " ttgctgacag\n",
            "      721 ctactatatc actactccat ctagtagtgg ccacgcccta tgaggcatat"
            " cctatcggaa\n",
            "      781 aacaataccc cccagtggca agagtcaatg aatcgtttac atttcaaatt"
            " tccaatgata\n",
            "      841 cctataaatc gtctgtagac aagacagctc aaataacata caattgcttc"
            " gacttaccga\n",
            "      901 gctggctttc gtttgactct agttctagaa cgttctcagg tgaaccttct"
            " tctgacttac\n",
            "      961 tatctgatgc gaacaccacg ttgtatttca atgtaatact cgagggtacg"
            " gactctgccg\n",
            "     1021 acagcacgtc tttgaacaat acataccaat ttgttgttac aaaccgtcca"
            " tccatctcgc\n",
            "     1081 tatcgtcaga tttcaatcta ttggcgttgt taaaaaacta tggttatact"
            " aacggcaaaa\n",
            "     1141 acgctctgaa actagatcct aatgaagtct tcaacgtgac ttttgaccgt"
            " tcaatgttca\n",
            "     1201 ctaacgaaga atccattgtg tcgtattacg gacgttctca gttgtataat"
            " gcgccgttac\n",
            "     1261 ccaattggct gttcttcgat tctggcgagt tgaagtttac tgggacggca"
            " ccggtgataa\n",
            "     1321 actcggcgat tgctccagaa acaagctaca gttttgtcat catcgctaca"
            " gacattgaag\n",
            "     1381 gattttctgc cgttgaggta gaattcgaat tagtcatcgg ggctcaccag"
            " ttaactacct\n",
            "     1441 ctattcaaaa tagtttgata atcaacgtta ctgacacagg taacgtttca"
            " tatgacttac\n",
            "     1501 ctctaaacta tgtttatctc gatgacgatc ctatttcttc tgataaattg"
            " ggttctataa\n",
            "     1561 acttattgga tgctccagac tgggtggcat tagataatgc taccatttcc"
            " gggtctgtcc\n",
            "     1621 cagatgaatt actcggtaag aactccaatc ctgccaattt ttctgtgtcc"
            " atttatgata\n",
            "     1681 cttatggtga tgtgatttat ttcaacttcg aagttgtctc cacaacggat"
            " ttgtttgcca\n",
            "     1741 ttagttctct tcccaatatt aacgctacaa ggggtgaatg gttctcctac"
            " tattttttgc\n",
            "     1801 cttctcagtt tacagactac gtgaatacaa acgtttcatt agagtttact"
            " aattcaagcc\n",
            "     1861 aagaccatga ctgggtgaaa ttccaatcat ctaatttaac attagctgga"
            " gaagtgccca\n",
            "     1921 agaatttcga caagctttca ttaggtttga aagcgaacca aggttcacaa"
            " tctcaagagc\n",
            "     1981 tatattttaa catcattggc atggattcaa agataactca ctcaaaccac"
            " agtgcgaatg\n",
            "     2041 caacgtccac aagaagttct caccactcca cctcaacaag ttcttacaca"
            " tcttctactt\n",
            "     2101 acactgcaaa aatttcttct acctccgctg ctgctacttc ttctgctcca"
            " gcagcgctgc\n",
            "     2161 cagcagccaa taaaacttca tctcacaata aaaaagcagt agcaattgcg"
            " tgcggtgttg\n",
            "     2221 ctatcccatt aggcgttatc ctagtagctc tcatttgctt cctaatattc"
            " tggagacgca\n",
            "     2281 gaagggaaaa tccagacgat gaaaacttac cgcatgctat tagtggacct"
            " gatttgaata\n",
            "     2341 atcctgcaaa taaaccaaat caagaaaacg ctacaccttt gaacaacccc"
            " tttgatgatg\n",
            "     2401 atgcttcctc gtacgatgat acttcaatag caagaagatt ggctgctttg"
            " aacactttga\n",
            "     2461 aattggataa ccactctgcc actgaatctg atatttccag cgtggatgaa"
            " aagagagatt\n",
            "     2521 ctctatcagg tatgaataca tacaatgatc agttccaatc ccaaagtaaa"
            " gaagaattat\n",
            "     2581 tagcaaaacc cccagtacag cctccagaga gcccgttctt tgacccacag"
            " aataggtctt\n",
            "     2641 cttctgtgta tatggatagt gaaccagcag taaataaatc ctggcgatat"
            " actggcaacc\n",
            "     2701 tgtcaccagt ctctgatatt gtcagagaca gttacggatc acaaaaaact"
            " gttgatacag\n",
            "     2761 aaaaactttt cgatttagaa gcaccagaga aggaaaaacg tacgtcaagg"
            " gatgtcacta\n",
            "     2821 tgtcttcact ggacccttgg aacagcaata ttagcccttc tcccgtaaga"
            " aaatcagtaa\n",
            "     2881 caccatcacc atataacgta acgaagcatc gtaaccgcca cttacaaaat"
            " attcaagact\n",
            "     2941 ctcaaagcgg taaaaacgga atcactccca caacaatgtc aacttcatct"
            " tctgacgatt\n",
            "     3001 ttgttccggt taaagatggt gaaaattttt gctgggtcca tagcatggaa"
            " ccagacagaa\n",
            "     3061 gaccaagtaa gaaaaggtta gtagattttt caaataagag taatgtcaat"
            " gttggtcaag\n",
            "     3121 ttaaggacat tcacggacgc atcccagaaa tgctgtgatt atacgcaacg"
            " atattttgct\n",
            "     3181 taattttatt ttcctgtttt attttttatt agtggtttac agatacccta"
            " tattttattt\n",
            "     3241 agtttttata cttagagaca tttaatttta attccattct tcaaatttca"
            " tttttgcact\n",
            "     3301 taaaacaaag atccaaaaat gctctcgccc tcttcatatt gagaatacac"
            " tccattcaaa\n",
            "     3361 attttgtcgt caccgctgat taatttttca ctaaactgat gaataatcaa"
            " aggccccacg\n",
            "     3421 tcagaaccga ctaaagaagt gagttttatt ttaggaggtt gaaaaccatt"
            " attgtctggt\n",
            "     3481 aaattttcat cttcttgaca tttaacccag tttgaatccc tttcaatttc"
            " tgctttttcc\n",
            "     3541 tccaaactat cgaccctcct gtttctgtcc aacttatgtc ctagttccaa"
            " ttcgatcgca\n",
            "     3601 ttaataactg cttcaaatgt tattgtgtca tcgttgactt taggtaattt"
            " ctccaaatgc\n",
            "     3661 ataatcaaac tatttaagga agatcggaat tcgtcgaaca cttcagtttc"
            " cgtaatgatc\n",
            "     3721 tgatcgtctt tatccacatg ttgtaattca ctaaaatcta aaacgtattt"
            " ttcaatgcat\n",
            "     3781 aaatcgttct ttttattaat aatgcagatg gaaaatctgt aaacgtgcgt"
            " taatttagaa\n",
            "     3841 agaacatcca gtataagttc ttctatatag tcaattaaag caggatgcct"
            " attaatggga\n",
            "     3901 acgaactgcg gcaagttgaa tgactggtaa gtagtgtagt cgaatgactg"
            " aggtgggtat\n",
            "     3961 acatttctat aaaataaaat caaattaatg tagcatttta agtataccct"
            " cagccacttc\n",
            "     4021 tctacccatc tattcataaa gctgacgcaa cgattactat tttttttttc"
            " ttcttggatc\n",
            "     4081 tcagtcgtcg caaaaacgta taccttcttt ttccgacctt ttttttagct"
            " ttctggaaaa\n",
            "     4141 gtttatatta gttaaacagg gtctagtctt agtgtgaaag ctagtggttt"
            " cgattgactg\n",
            "     4201 atattaagaa agtggaaatt aaattagtag tgtagacgta tatgcatatg"
            " tatttctcgc\n",
            "     4261 ctgtttatgt ttctacgtac ttttgattta tagcaagggg aaaagaaata"
            " catactattt\n",
            "     4321 tttggtaaag gtgaaagcat aatgtaaaag ctagaataaa atggacgaaa"
            " taaagagagg\n",
            "     4381 cttagttcat cttttttcca aaaagcaccc aatgataata actaaaatga"
            " aaaggatttg\n",
            "     4441 ccatctgtca gcaacatcag ttgtgtgagc aataataaaa tcatcacctc"
            " cgttgccttt\n",
            "     4501 agcgcgtttg tcgtttgtat cttccgtaat tttagtctta tcaatgggaa"
            " tcataaattt\n",
            "     4561 tccaatgaat tagcaatttc gtccaattct ttttgagctt cttcatattt"
            " gctttggaat\n",
            "     4621 tcttcgcact tcttttccca ttcatctctt tcttcttcca aagcaacgat"
            " ccttctaccc\n",
            "     4681 atttgctcag agttcaaatc ggcctctttc agtttatcca ttgcttcctt"
            " cagtttggct\n",
            "     4741 tcactgtctt ctagctgttg ttctagatcc tggtttttct tggtgtagtt"
            " ctcattatta\n",
            "     4801 gatctcaagt tattggagtc ttcagccaat tgctttgtat cagacaattg"
            " actctctaac\n",
            "     4861 ttctccactt cactgtcgag ttgctcgttt ttagcggaca aagatttaat"
            " ctcgttttct\n",
            "     4921 ttttcagtgt tagattgctc taattctttg agctgttctc tcagctcctc"
            " atatttttct\n",
            "     4981 tgccatgact cagattctaa ttttaagcta ttcaatttct ctttgatc"
            "\n",
            "//\n"
        ]}
        for ext in ('fna', 'gbk'):
            output_gm_fp = join(self.working_dir, "id." + ext)
            with open(output_gm_fp, 'r') as output_gm_f:
                reformat_gm_act = output_gm_f.readlines()
            self.assertListEqual(reformat_gm_exp[ext], reformat_gm_act)

    def test_reformat_egid(self):
        """ Test functionality of reformat_egid()
        """
        reformat_egid(self.species_genome_fp,
                      self.working_dir)
        reformat_egid_exp = {'fna': [
            ">locus001\n",
            "GATCCTCCATATACAACGGTATCTCCACCTCAGGTTTAGATCTCAACAACGGAACCATTGCCGA"
            "CATGAGACAGTTAGGTATCGTCGAGAGTTACAAGCTAAAACGAGCAGTAGTCAGCTCTGCATCT"
            "GAAGCCGCTGAAGTTCTACTAAGGGTGGATAACATCATCCGTGCAAGACCAAGAACCGCCAATA"
            "GACAACATATGTAACATATTTAGGATATACCTCGAAAATAATAAACCGCCACACTGTCATTATT"
            "ATAATTAGAAACAGAACGCAAAAATTATCCACTATATAATTCAAAGACGCGAAAAAAAAAGAAC"
            "AACGCGTCATAGAACTTTTGGCAATTCGCGTCACAAATAAATTTTGGCAACTTATGTTTCCTCT"
            "TCGAGCAGTACTCGAGCCCTGTCTCAAGAATGTAATAATACCCATCGTAGGTATGGTTAAAGAT"
            "AGCATCTCCACAACCTCAAAGCTCCTTGCCGAGAGTCGCCCTCCTTTGTCGAGTAATTTTCACT"
            "TTTCATATGAGAACTTATTTTCTTATTCTTTACTCTCACATCCTGTAGTGATTGACACTGCAAC"
            "AGCCACCATCACTAGAAGAACAGAACAATTACTTAATAGAAAAATTATATCTTCCTCGAAACGA"
            "TTTCCTGCTTCCAACATCTACGTATATCAAGAAGCATTCACTTACCATGACACAGCTTCAGATT"
            "TCATTATTGCTGACAGCTACTATATCACTACTCCATCTAGTAGTGGCCACGCCCTATGAGGCAT"
            "ATCCTATCGGAAAACAATACCCCCCAGTGGCAAGAGTCAATGAATCGTTTACATTTCAAATTTC"
            "CAATGATACCTATAAATCGTCTGTAGACAAGACAGCTCAAATAACATACAATTGCTTCGACTTA"
            "CCGAGCTGGCTTTCGTTTGACTCTAGTTCTAGAACGTTCTCAGGTGAACCTTCTTCTGACTTAC"
            "TATCTGATGCGAACACCACGTTGTATTTCAATGTAATACTCGAGGGTACGGACTCTGCCGACAG"
            "CACGTCTTTGAACAATACATACCAATTTGTTGTTACAAACCGTCCATCCATCTCGCTATCGTCA"
            "GATTTCAATCTATTGGCGTTGTTAAAAAACTATGGTTATACTAACGGCAAAAACGCTCTGAAAC"
            "TAGATCCTAATGAAGTCTTCAACGTGACTTTTGACCGTTCAATGTTCACTAACGAAGAATCCAT"
            "TGTGTCGTATTACGGACGTTCTCAGTTGTATAATGCGCCGTTACCCAATTGGCTGTTCTTCGAT"
            "TCTGGCGAGTTGAAGTTTACTGGGACGGCACCGGTGATAAACTCGGCGATTGCTCCAGAAACAA"
            "GCTACAGTTTTGTCATCATCGCTACAGACATTGAAGGATTTTCTGCCGTTGAGGTAGAATTCGA"
            "ATTAGTCATCGGGGCTCACCAGTTAACTACCTCTATTCAAAATAGTTTGATAATCAACGTTACT"
            "GACACAGGTAACGTTTCATATGACTTACCTCTAAACTATGTTTATCTCGATGACGATCCTATTT"
            "CTTCTGATAAATTGGGTTCTATAAACTTATTGGATGCTCCAGACTGGGTGGCATTAGATAATGC"
            "TACCATTTCCGGGTCTGTCCCAGATGAATTACTCGGTAAGAACTCCAATCCTGCCAATTTTTCT"
            "GTGTCCATTTATGATACTTATGGTGATGTGATTTATTTCAACTTCGAAGTTGTCTCCACAACGG"
            "ATTTGTTTGCCATTAGTTCTCTTCCCAATATTAACGCTACAAGGGGTGAATGGTTCTCCTACTA"
            "TTTTTTGCCTTCTCAGTTTACAGACTACGTGAATACAAACGTTTCATTAGAGTTTACTAATTCA"
            "AGCCAAGACCATGACTGGGTGAAATTCCAATCATCTAATTTAACATTAGCTGGAGAAGTGCCCA"
            "AGAATTTCGACAAGCTTTCATTAGGTTTGAAAGCGAACCAAGGTTCACAATCTCAAGAGCTATA"
            "TTTTAACATCATTGGCATGGATTCAAAGATAACTCACTCAAACCACAGTGCGAATGCAACGTCC"
            "ACAAGAAGTTCTCACCACTCCACCTCAACAAGTTCTTACACATCTTCTACTTACACTGCAAAAA"
            "TTTCTTCTACCTCCGCTGCTGCTACTTCTTCTGCTCCAGCAGCGCTGCCAGCAGCCAATAAAAC"
            "TTCATCTCACAATAAAAAAGCAGTAGCAATTGCGTGCGGTGTTGCTATCCCATTAGGCGTTATC"
            "CTAGTAGCTCTCATTTGCTTCCTAATATTCTGGAGACGCAGAAGGGAAAATCCAGACGATGAAA"
            "ACTTACCGCATGCTATTAGTGGACCTGATTTGAATAATCCTGCAAATAAACCAAATCAAGAAAA"
            "CGCTACACCTTTGAACAACCCCTTTGATGATGATGCTTCCTCGTACGATGATACTTCAATAGCA"
            "AGAAGATTGGCTGCTTTGAACACTTTGAAATTGGATAACCACTCTGCCACTGAATCTGATATTT"
            "CCAGCGTGGATGAAAAGAGAGATTCTCTATCAGGTATGAATACATACAATGATCAGTTCCAATC"
            "CCAAAGTAAAGAAGAATTATTAGCAAAACCCCCAGTACAGCCTCCAGAGAGCCCGTTCTTTGAC"
            "CCACAGAATAGGTCTTCTTCTGTGTATATGGATAGTGAACCAGCAGTAAATAAATCCTGGCGAT"
            "ATACTGGCAACCTGTCACCAGTCTCTGATATTGTCAGAGACAGTTACGGATCACAAAAAACTGT"
            "TGATACAGAAAAACTTTTCGATTTAGAAGCACCAGAGAAGGAAAAACGTACGTCAAGGGATGTC"
            "ACTATGTCTTCACTGGACCCTTGGAACAGCAATATTAGCCCTTCTCCCGTAAGAAAATCAGTAA"
            "CACCATCACCATATAACGTAACGAAGCATCGTAACCGCCACTTACAAAATATTCAAGACTCTCA"
            "AAGCGGTAAAAACGGAATCACTCCCACAACAATGTCAACTTCATCTTCTGACGATTTTGTTCCG"
            "GTTAAAGATGGTGAAAATTTTTGCTGGGTCCATAGCATGGAACCAGACAGAAGACCAAGTAAGA"
            "AAAGGTTAGTAGATTTTTCAAATAAGAGTAATGTCAATGTTGGTCAAGTTAAGGACATTCACGG"
            "ACGCATCCCAGAAATGCTGTGATTATACGCAACGATATTTTGCTTAATTTTATTTTCCTGTTTT"
            "ATTTTTTATTAGTGGTTTACAGATACCCTATATTTTATTTAGTTTTTATACTTAGAGACATTTA"
            "ATTTTAATTCCATTCTTCAAATTTCATTTTTGCACTTAAAACAAAGATCCAAAAATGCTCTCGC"
            "CCTCTTCATATTGAGAATACACTCCATTCAAAATTTTGTCGTCACCGCTGATTAATTTTTCACT"
            "AAACTGATGAATAATCAAAGGCCCCACGTCAGAACCGACTAAAGAAGTGAGTTTTATTTTAGGA"
            "GGTTGAAAACCATTATTGTCTGGTAAATTTTCATCTTCTTGACATTTAACCCAGTTTGAATCCC"
            "TTTCAATTTCTGCTTTTTCCTCCAAACTATCGACCCTCCTGTTTCTGTCCAACTTATGTCCTAG"
            "TTCCAATTCGATCGCATTAATAACTGCTTCAAATGTTATTGTGTCATCGTTGACTTTAGGTAAT"
            "TTCTCCAAATGCATAATCAAACTATTTAAGGAAGATCGGAATTCGTCGAACACTTCAGTTTCCG"
            "TAATGATCTGATCGTCTTTATCCACATGTTGTAATTCACTAAAATCTAAAACGTATTTTTCAAT"
            "GCATAAATCGTTCTTTTTATTAATAATGCAGATGGAAAATCTGTAAACGTGCGTTAATTTAGAA"
            "AGAACATCCAGTATAAGTTCTTCTATATAGTCAATTAAAGCAGGATGCCTATTAATGGGAACGA"
            "ACTGCGGCAAGTTGAATGACTGGTAAGTAGTGTAGTCGAATGACTGAGGTGGGTATACATTTCT"
            "ATAAAATAAAATCAAATTAATGTAGCATTTTAAGTATACCCTCAGCCACTTCTCTACCCATCTA"
            "TTCATAAAGCTGACGCAACGATTACTATTTTTTTTTTCTTCTTGGATCTCAGTCGTCGCAAAAA"
            "CGTATACCTTCTTTTTCCGACCTTTTTTTTAGCTTTCTGGAAAAGTTTATATTAGTTAAACAGG"
            "GTCTAGTCTTAGTGTGAAAGCTAGTGGTTTCGATTGACTGATATTAAGAAAGTGGAAATTAAAT"
            "TAGTAGTGTAGACGTATATGCATATGTATTTCTCGCCTGTTTATGTTTCTACGTACTTTTGATT"
            "TATAGCAAGGGGAAAAGAAATACATACTATTTTTTGGTAAAGGTGAAAGCATAATGTAAAAGCT"
            "AGAATAAAATGGACGAAATAAAGAGAGGCTTAGTTCATCTTTTTTCCAAAAAGCACCCAATGAT"
            "AATAACTAAAATGAAAAGGATTTGCCATCTGTCAGCAACATCAGTTGTGTGAGCAATAATAAAA"
            "TCATCACCTCCGTTGCCTTTAGCGCGTTTGTCGTTTGTATCTTCCGTAATTTTAGTCTTATCAA"
            "TGGGAATCATAAATTTTCCAATGAATTAGCAATTTCGTCCAATTCTTTTTGAGCTTCTTCATAT"
            "TTGCTTTGGAATTCTTCGCACTTCTTTTCCCATTCATCTCTTTCTTCTTCCAAAGCAACGATCC"
            "TTCTACCCATTTGCTCAGAGTTCAAATCGGCCTCTTTCAGTTTATCCATTGCTTCCTTCAGTTT"
            "GGCTTCACTGTCTTCTAGCTGTTGTTCTAGATCCTGGTTTTTCTTGGTGTAGTTCTCATTATTA"
            "GATCTCAAGTTATTGGAGTCTTCAGCCAATTGCTTTGTATCAGACAATTGACTCTCTAACTTCT"
            "CCACTTCACTGTCGAGTTGCTCGTTTTTAGCGGACAAAGATTTAATCTCGTTTTCTTTTTCAGT"
            "GTTAGATTGCTCTAATTCTTTGAGCTGTTCTCTCAGCTCCTCATATTTTTCTTGCCATGACTCA"
            "GATTCTAATTTTAAGCTATTCAATTTCTCTTTGATC\n"
        ], 'faa': [
            ">AAA98665.1\n",
            "SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEAAEVLLRVDNIIRARPRTANR"
            "QHM\n",
            ">AAA98666.1\n",
            "MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESFTFQISNDTYKSSVDKTAQIT"
            "YNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLYFNVILEGTDSADSTSLNNTYQFVVTNRP"
            "SISLSSDFNLLALLKNYGYTNGKNALKLDPNEVFNVTFDRSMFTNEESIVSYYGRSQLYNAPLP"
            "NWLFFDSGELKFTGTAPVINSAIAPETSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNS"
            "LIINVTDTGNVSYDLPLNYVYLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNS"
            "NPANFSVSIYDTYGDVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVS"
            "LEFTNSSQDHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNH"
            "SANATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIACGVA"
            "IPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLNNPFDDDASSY"
            "DDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQSQSKEELLAKPPVQPP"
            "ESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDSYGSQKTVDTEKLFDLEAPEKEK"
            "RTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSS"
            "SDDFVPVKDGENFCWVHSMEPDRRPSKKRLVDFSNKSNVNVGQVKDIHGRIPEML\n",
            ">AAA98667.1\n",
            "MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQFVPINRHPALIDYIEELILD"
            "VLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVDKDDQIITETEVFDEFRSSLNSLIMHL"
            "EKLPKVNDDTITFEAVINAIELELGHKLDRNRRVDSLEEKAEIERDSNWVKCQEDENLPDNNGF"
            "QPPKIKLTSLVGSDVGPLIIHQFSEKLISGDDKILNGVYSQYEEGESIFGSLF\n"
        ], 'ffn': [
            ">locus001:1-206\n",
            "GATCCTCCATATACAACGGTATCTCCACCTCAGGTTTAGATCTCAACAACGGAACCATTGCCGA"
            "CATGAGACAGTTAGGTATCGTCGAGAGTTACAAGCTAAAACGAGCAGTAGTCAGCTCTGCATCT"
            "GAAGCCGCTGAAGTTCTACTAAGGGTGGATAACATCATCCGTGCAAGACCAAGAACCGCCAATA"
            "GACAACATATGTAA\n",
            ">locus001:687-3158\n",
            "ATGACACAGCTTCAGATTTCATTATTGCTGACAGCTACTATATCACTACTCCATCTAGTAGTGG"
            "CCACGCCCTATGAGGCATATCCTATCGGAAAACAATACCCCCCAGTGGCAAGAGTCAATGAATC"
            "GTTTACATTTCAAATTTCCAATGATACCTATAAATCGTCTGTAGACAAGACAGCTCAAATAACA"
            "TACAATTGCTTCGACTTACCGAGCTGGCTTTCGTTTGACTCTAGTTCTAGAACGTTCTCAGGTG"
            "AACCTTCTTCTGACTTACTATCTGATGCGAACACCACGTTGTATTTCAATGTAATACTCGAGGG"
            "TACGGACTCTGCCGACAGCACGTCTTTGAACAATACATACCAATTTGTTGTTACAAACCGTCCA"
            "TCCATCTCGCTATCGTCAGATTTCAATCTATTGGCGTTGTTAAAAAACTATGGTTATACTAACG"
            "GCAAAAACGCTCTGAAACTAGATCCTAATGAAGTCTTCAACGTGACTTTTGACCGTTCAATGTT"
            "CACTAACGAAGAATCCATTGTGTCGTATTACGGACGTTCTCAGTTGTATAATGCGCCGTTACCC"
            "AATTGGCTGTTCTTCGATTCTGGCGAGTTGAAGTTTACTGGGACGGCACCGGTGATAAACTCGG"
            "CGATTGCTCCAGAAACAAGCTACAGTTTTGTCATCATCGCTACAGACATTGAAGGATTTTCTGC"
            "CGTTGAGGTAGAATTCGAATTAGTCATCGGGGCTCACCAGTTAACTACCTCTATTCAAAATAGT"
            "TTGATAATCAACGTTACTGACACAGGTAACGTTTCATATGACTTACCTCTAAACTATGTTTATC"
            "TCGATGACGATCCTATTTCTTCTGATAAATTGGGTTCTATAAACTTATTGGATGCTCCAGACTG"
            "GGTGGCATTAGATAATGCTACCATTTCCGGGTCTGTCCCAGATGAATTACTCGGTAAGAACTCC"
            "AATCCTGCCAATTTTTCTGTGTCCATTTATGATACTTATGGTGATGTGATTTATTTCAACTTCG"
            "AAGTTGTCTCCACAACGGATTTGTTTGCCATTAGTTCTCTTCCCAATATTAACGCTACAAGGGG"
            "TGAATGGTTCTCCTACTATTTTTTGCCTTCTCAGTTTACAGACTACGTGAATACAAACGTTTCA"
            "TTAGAGTTTACTAATTCAAGCCAAGACCATGACTGGGTGAAATTCCAATCATCTAATTTAACAT"
            "TAGCTGGAGAAGTGCCCAAGAATTTCGACAAGCTTTCATTAGGTTTGAAAGCGAACCAAGGTTC"
            "ACAATCTCAAGAGCTATATTTTAACATCATTGGCATGGATTCAAAGATAACTCACTCAAACCAC"
            "AGTGCGAATGCAACGTCCACAAGAAGTTCTCACCACTCCACCTCAACAAGTTCTTACACATCTT"
            "CTACTTACACTGCAAAAATTTCTTCTACCTCCGCTGCTGCTACTTCTTCTGCTCCAGCAGCGCT"
            "GCCAGCAGCCAATAAAACTTCATCTCACAATAAAAAAGCAGTAGCAATTGCGTGCGGTGTTGCT"
            "ATCCCATTAGGCGTTATCCTAGTAGCTCTCATTTGCTTCCTAATATTCTGGAGACGCAGAAGGG"
            "AAAATCCAGACGATGAAAACTTACCGCATGCTATTAGTGGACCTGATTTGAATAATCCTGCAAA"
            "TAAACCAAATCAAGAAAACGCTACACCTTTGAACAACCCCTTTGATGATGATGCTTCCTCGTAC"
            "GATGATACTTCAATAGCAAGAAGATTGGCTGCTTTGAACACTTTGAAATTGGATAACCACTCTG"
            "CCACTGAATCTGATATTTCCAGCGTGGATGAAAAGAGAGATTCTCTATCAGGTATGAATACATA"
            "CAATGATCAGTTCCAATCCCAAAGTAAAGAAGAATTATTAGCAAAACCCCCAGTACAGCCTCCA"
            "GAGAGCCCGTTCTTTGACCCACAGAATAGGTCTTCTTCTGTGTATATGGATAGTGAACCAGCAG"
            "TAAATAAATCCTGGCGATATACTGGCAACCTGTCACCAGTCTCTGATATTGTCAGAGACAGTTA"
            "CGGATCACAAAAAACTGTTGATACAGAAAAACTTTTCGATTTAGAAGCACCAGAGAAGGAAAAA"
            "CGTACGTCAAGGGATGTCACTATGTCTTCACTGGACCCTTGGAACAGCAATATTAGCCCTTCTC"
            "CCGTAAGAAAATCAGTAACACCATCACCATATAACGTAACGAAGCATCGTAACCGCCACTTACA"
            "AAATATTCAAGACTCTCAAAGCGGTAAAAACGGAATCACTCCCACAACAATGTCAACTTCATCT"
            "TCTGACGATTTTGTTCCGGTTAAAGATGGTGAAAATTTTTGCTGGGTCCATAGCATGGAACCAG"
            "ACAGAAGACCAAGTAAGAAAAGGTTAGTAGATTTTTCAAATAAGAGTAATGTCAATGTTGGTCA"
            "AGTTAAGGACATTCACGGACGCATCCCAGAAATGCTGTGA\n",
            ">locus001:c4037-3300\n",
            "ATGAATAGATGGGTAGAGAAGTGGCTGAGGGTATACTTAAAATGCTACATTAATTTGATTTTAT"
            "TTTATAGAAATGTATACCCACCTCAGTCATTCGACTACACTACTTACCAGTCATTCAACTTGCC"
            "GCAGTTCGTTCCCATTAATAGGCATCCTGCTTTAATTGACTATATAGAAGAACTTATACTGGAT"
            "GTTCTTTCTAAATTAACGCACGTTTACAGATTTTCCATCTGCATTATTAATAAAAAGAACGATT"
            "TATGCATTGAAAAATACGTTTTAGATTTTAGTGAATTACAACATGTGGATAAAGACGATCAGAT"
            "CATTACGGAAACTGAAGTGTTCGACGAATTCCGATCTTCCTTAAATAGTTTGATTATGCATTTG"
            "GAGAAATTACCTAAAGTCAACGATGACACAATAACATTTGAAGCAGTTATTAATGCGATCGAAT"
            "TGGAACTAGGACATAAGTTGGACAGAAACAGGAGGGTCGATAGTTTGGAGGAAAAAGCAGAAAT"
            "TGAAAGGGATTCAAACTGGGTTAAATGTCAAGAAGATGAAAATTTACCAGACAATAATGGTTTT"
            "CAACCTCCTAAAATAAAACTCACTTCTTTAGTCGGTTCTGACGTGGGGCCTTTGATTATTCATC"
            "AGTTTAGTGAAAAATTAATCAGCGGTGACGACAAAATTTTGAATGGAGTGTATTCTCAATATGA"
            "AGAGGGCGAGAGCATTTTTGGATCTTTGTTTTAA\n"
        ], 'ptt': [
            "locus001\n",
            "3 proteins\n",
            "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct"
            "\n",
            "1..206	+	67	1	-	gene1	-	-	-\n",
            "687..3158	+	823	2	-	gene2	-	-	-\n",
            "3300..4037	-	245	3	-	gene3	-	-	-\n"
        ], 'gbk': [
            "LOCUS       locus001   5028 bp   DNA   circular   CON   01-JAN-1"
            "900\n",
            "FEATURES             Location/Qualifiers\n",
            "     gene            1..206\n",
            "                     /locus_tag=gene1\n",
            "     CDS             1..206\n",
            "                     /locus_tag=gene1\n",
            "                     /protein_id=AAA98665.1\n",
            "                     /translation=SSIYNGISTSGLDLNNGTIADMRQLGIVES"
            "YKLKRAVVSSASEAAEVLLRVDNIIRARPRTANRQHM\n",
            "     gene            687..3158\n",
            "                     /locus_tag=gene2\n",
            "     CDS             687..3158\n",
            "                     /locus_tag=gene2\n",
            "                     /protein_id=AAA98666.1\n",
            "                     /translation=MTQLQISLLLTATISLLHLVVATPYEAYPI"
            "GKQYPPVARVNESFTFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSD"
            "ANTTLYFNVILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDP"
            "NEVFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPETSYS"
            "FVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYVYLDDDPISSD"
            "KLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYGDVIYFNFEVVSTTDLF"
            "AISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQDHDWVKFQSSNLTLAGEVPKNF"
            "DKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSANATSTRSSHHSTSTSSYTSSTYTAKISS"
            "TSAAATSSAPAALPAANKTSSHNKKAVAIACGVAIPLGVILVALICFLIFWRRRRENPDDENLP"
            "HAISGPDLNNPANKPNQENATPLNNPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSV"
            "DEKRDSLSGMNTYNDQFQSQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTG"
            "NLSPVSDIVRDSYGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPS"
            "PYNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRL"
            "VDFSNKSNVNVGQVKDIHGRIPEML\n",
            "     gene            complement(3300..4037)\n",
            "                     /locus_tag=gene3\n",
            "     CDS             complement(3300..4037)\n",
            "                     /locus_tag=gene3\n",
            "                     /protein_id=AAA98667.1\n",
            "                     /translation=MNRWVEKWLRVYLKCYINLILFYRNVYPPQ"
            "SFDYTTYQSFNLPQFVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLD"
            "FSELQHVDKDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDR"
            "NRRVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEKLISG"
            "DDKILNGVYSQYEEGESIFGSLF\n",
            "ORIGIN\n",
            "        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac"
            " ggaaccattg\n",
            "       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca"
            " gtagtcagct\n",
            "      121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt"
            " gcaagaccaa\n",
            "      181 gaaccgccaa tagacaacat atgtaacata tttaggatat acctcgaaaa"
            " taataaaccg\n",
            "      241 ccacactgtc attattataa ttagaaacag aacgcaaaaa ttatccacta"
            " tataattcaa\n",
            "      301 agacgcgaaa aaaaaagaac aacgcgtcat agaacttttg gcaattcgcg"
            " tcacaaataa\n",
            "      361 attttggcaa cttatgtttc ctcttcgagc agtactcgag ccctgtctca"
            " agaatgtaat\n",
            "      421 aatacccatc gtaggtatgg ttaaagatag catctccaca acctcaaagc"
            " tccttgccga\n",
            "      481 gagtcgccct cctttgtcga gtaattttca cttttcatat gagaacttat"
            " tttcttattc\n",
            "      541 tttactctca catcctgtag tgattgacac tgcaacagcc accatcacta"
            " gaagaacaga\n",
            "      601 acaattactt aatagaaaaa ttatatcttc ctcgaaacga tttcctgctt"
            " ccaacatcta\n",
            "      661 cgtatatcaa gaagcattca cttaccatga cacagcttca gatttcatta"
            " ttgctgacag\n",
            "      721 ctactatatc actactccat ctagtagtgg ccacgcccta tgaggcatat"
            " cctatcggaa\n",
            "      781 aacaataccc cccagtggca agagtcaatg aatcgtttac atttcaaatt"
            " tccaatgata\n",
            "      841 cctataaatc gtctgtagac aagacagctc aaataacata caattgcttc"
            " gacttaccga\n",
            "      901 gctggctttc gtttgactct agttctagaa cgttctcagg tgaaccttct"
            " tctgacttac\n",
            "      961 tatctgatgc gaacaccacg ttgtatttca atgtaatact cgagggtacg"
            " gactctgccg\n",
            "     1021 acagcacgtc tttgaacaat acataccaat ttgttgttac aaaccgtcca"
            " tccatctcgc\n",
            "     1081 tatcgtcaga tttcaatcta ttggcgttgt taaaaaacta tggttatact"
            " aacggcaaaa\n",
            "     1141 acgctctgaa actagatcct aatgaagtct tcaacgtgac ttttgaccgt"
            " tcaatgttca\n",
            "     1201 ctaacgaaga atccattgtg tcgtattacg gacgttctca gttgtataat"
            " gcgccgttac\n",
            "     1261 ccaattggct gttcttcgat tctggcgagt tgaagtttac tgggacggca"
            " ccggtgataa\n",
            "     1321 actcggcgat tgctccagaa acaagctaca gttttgtcat catcgctaca"
            " gacattgaag\n",
            "     1381 gattttctgc cgttgaggta gaattcgaat tagtcatcgg ggctcaccag"
            " ttaactacct\n",
            "     1441 ctattcaaaa tagtttgata atcaacgtta ctgacacagg taacgtttca"
            " tatgacttac\n",
            "     1501 ctctaaacta tgtttatctc gatgacgatc ctatttcttc tgataaattg"
            " ggttctataa\n",
            "     1561 acttattgga tgctccagac tgggtggcat tagataatgc taccatttcc"
            " gggtctgtcc\n",
            "     1621 cagatgaatt actcggtaag aactccaatc ctgccaattt ttctgtgtcc"
            " atttatgata\n",
            "     1681 cttatggtga tgtgatttat ttcaacttcg aagttgtctc cacaacggat"
            " ttgtttgcca\n",
            "     1741 ttagttctct tcccaatatt aacgctacaa ggggtgaatg gttctcctac"
            " tattttttgc\n",
            "     1801 cttctcagtt tacagactac gtgaatacaa acgtttcatt agagtttact"
            " aattcaagcc\n",
            "     1861 aagaccatga ctgggtgaaa ttccaatcat ctaatttaac attagctgga"
            " gaagtgccca\n",
            "     1921 agaatttcga caagctttca ttaggtttga aagcgaacca aggttcacaa"
            " tctcaagagc\n",
            "     1981 tatattttaa catcattggc atggattcaa agataactca ctcaaaccac"
            " agtgcgaatg\n",
            "     2041 caacgtccac aagaagttct caccactcca cctcaacaag ttcttacaca"
            " tcttctactt\n",
            "     2101 acactgcaaa aatttcttct acctccgctg ctgctacttc ttctgctcca"
            " gcagcgctgc\n",
            "     2161 cagcagccaa taaaacttca tctcacaata aaaaagcagt agcaattgcg"
            " tgcggtgttg\n",
            "     2221 ctatcccatt aggcgttatc ctagtagctc tcatttgctt cctaatattc"
            " tggagacgca\n",
            "     2281 gaagggaaaa tccagacgat gaaaacttac cgcatgctat tagtggacct"
            " gatttgaata\n",
            "     2341 atcctgcaaa taaaccaaat caagaaaacg ctacaccttt gaacaacccc"
            " tttgatgatg\n",
            "     2401 atgcttcctc gtacgatgat acttcaatag caagaagatt ggctgctttg"
            " aacactttga\n",
            "     2461 aattggataa ccactctgcc actgaatctg atatttccag cgtggatgaa"
            " aagagagatt\n",
            "     2521 ctctatcagg tatgaataca tacaatgatc agttccaatc ccaaagtaaa"
            " gaagaattat\n",
            "     2581 tagcaaaacc cccagtacag cctccagaga gcccgttctt tgacccacag"
            " aataggtctt\n",
            "     2641 cttctgtgta tatggatagt gaaccagcag taaataaatc ctggcgatat"
            " actggcaacc\n",
            "     2701 tgtcaccagt ctctgatatt gtcagagaca gttacggatc acaaaaaact"
            " gttgatacag\n",
            "     2761 aaaaactttt cgatttagaa gcaccagaga aggaaaaacg tacgtcaagg"
            " gatgtcacta\n",
            "     2821 tgtcttcact ggacccttgg aacagcaata ttagcccttc tcccgtaaga"
            " aaatcagtaa\n",
            "     2881 caccatcacc atataacgta acgaagcatc gtaaccgcca cttacaaaat"
            " attcaagact\n",
            "     2941 ctcaaagcgg taaaaacgga atcactccca caacaatgtc aacttcatct"
            " tctgacgatt\n",
            "     3001 ttgttccggt taaagatggt gaaaattttt gctgggtcca tagcatggaa"
            " ccagacagaa\n",
            "     3061 gaccaagtaa gaaaaggtta gtagattttt caaataagag taatgtcaat"
            " gttggtcaag\n",
            "     3121 ttaaggacat tcacggacgc atcccagaaa tgctgtgatt atacgcaacg"
            " atattttgct\n",
            "     3181 taattttatt ttcctgtttt attttttatt agtggtttac agatacccta"
            " tattttattt\n",
            "     3241 agtttttata cttagagaca tttaatttta attccattct tcaaatttca"
            " tttttgcact\n",
            "     3301 taaaacaaag atccaaaaat gctctcgccc tcttcatatt gagaatacac"
            " tccattcaaa\n",
            "     3361 attttgtcgt caccgctgat taatttttca ctaaactgat gaataatcaa"
            " aggccccacg\n",
            "     3421 tcagaaccga ctaaagaagt gagttttatt ttaggaggtt gaaaaccatt"
            " attgtctggt\n",
            "     3481 aaattttcat cttcttgaca tttaacccag tttgaatccc tttcaatttc"
            " tgctttttcc\n",
            "     3541 tccaaactat cgaccctcct gtttctgtcc aacttatgtc ctagttccaa"
            " ttcgatcgca\n",
            "     3601 ttaataactg cttcaaatgt tattgtgtca tcgttgactt taggtaattt"
            " ctccaaatgc\n",
            "     3661 ataatcaaac tatttaagga agatcggaat tcgtcgaaca cttcagtttc"
            " cgtaatgatc\n",
            "     3721 tgatcgtctt tatccacatg ttgtaattca ctaaaatcta aaacgtattt"
            " ttcaatgcat\n",
            "     3781 aaatcgttct ttttattaat aatgcagatg gaaaatctgt aaacgtgcgt"
            " taatttagaa\n",
            "     3841 agaacatcca gtataagttc ttctatatag tcaattaaag caggatgcct"
            " attaatggga\n",
            "     3901 acgaactgcg gcaagttgaa tgactggtaa gtagtgtagt cgaatgactg"
            " aggtgggtat\n",
            "     3961 acatttctat aaaataaaat caaattaatg tagcatttta agtataccct"
            " cagccacttc\n",
            "     4021 tctacccatc tattcataaa gctgacgcaa cgattactat tttttttttc"
            " ttcttggatc\n",
            "     4081 tcagtcgtcg caaaaacgta taccttcttt ttccgacctt ttttttagct"
            " ttctggaaaa\n",
            "     4141 gtttatatta gttaaacagg gtctagtctt agtgtgaaag ctagtggttt"
            " cgattgactg\n",
            "     4201 atattaagaa agtggaaatt aaattagtag tgtagacgta tatgcatatg"
            " tatttctcgc\n",
            "     4261 ctgtttatgt ttctacgtac ttttgattta tagcaagggg aaaagaaata"
            " catactattt\n",
            "     4321 tttggtaaag gtgaaagcat aatgtaaaag ctagaataaa atggacgaaa"
            " taaagagagg\n",
            "     4381 cttagttcat cttttttcca aaaagcaccc aatgataata actaaaatga"
            " aaaggatttg\n",
            "     4441 ccatctgtca gcaacatcag ttgtgtgagc aataataaaa tcatcacctc"
            " cgttgccttt\n",
            "     4501 agcgcgtttg tcgtttgtat cttccgtaat tttagtctta tcaatgggaa"
            " tcataaattt\n",
            "     4561 tccaatgaat tagcaatttc gtccaattct ttttgagctt cttcatattt"
            " gctttggaat\n",
            "     4621 tcttcgcact tcttttccca ttcatctctt tcttcttcca aagcaacgat"
            " ccttctaccc\n",
            "     4681 atttgctcag agttcaaatc ggcctctttc agtttatcca ttgcttcctt"
            " cagtttggct\n",
            "     4741 tcactgtctt ctagctgttg ttctagatcc tggtttttct tggtgtagtt"
            " ctcattatta\n",
            "     4801 gatctcaagt tattggagtc ttcagccaat tgctttgtat cagacaattg"
            " actctctaac\n",
            "     4861 ttctccactt cactgtcgag ttgctcgttt ttagcggaca aagatttaat"
            " ctcgttttct\n",
            "     4921 ttttcagtgt tagattgctc taattctttg agctgttctc tcagctcctc"
            " atatttttct\n",
            "     4981 tgccatgact cagattctaa ttttaagcta ttcaatttct ctttgatc"
            "\n",
            "//\n"
        ]}
        for ext in ('fna', 'faa', 'ffn', 'ptt', 'gbk'):
            output_egid_fp = join(self.working_dir, "id." + ext)
            with open(output_egid_fp, 'r') as output_egid_f:
                reformat_egid_act = output_egid_f.readlines()
            self.assertListEqual(reformat_egid_exp[ext], reformat_egid_act)


# 10 species
species_tree = ("(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
                "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
                "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
                "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
                "1.6606127,SE007:0.70000178):1.6331374):1.594016;")
# 10 species, 1 with same label
species_tree_2 = ("(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
                  "0.86305436):0.21024487,(SE001:0.56704221,SE009:0.5014676):"
                  "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
                  "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
                  "1.6606127,SE007:0.70000178):1.6331374):1.594016;")
# 10 species, 16 genes (gain)
gene_tree_1 = ("(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):"
               "3.7761166,SE008_01623:5.9256042):0.2102448,(SE006_01623:"
               "5.2329068,SE009_01623:5.2329068):0.9029422):0.2054233,"
               "SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):"
               "0.7293362,SE003_01623:7.4420648):0.2444784,((SE002_01623:"
               "6.0534057,SE007_01623:6.0534057):0.4589905,((((SE001_04123:"
               "2.1494876,SE010_04123:2.1494876):3.7761166,SE008_04123:"
               "5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_04123:"
               "5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):"
               "0.1711239):1.174147):1.594016;")
# 10 species, 10 genes
gene_tree_2 = ("(((((((SE001_00009:2.1494876,SE010_00009:2.1494876):"
               "3.7761166,SE008_00009:5.9256042):0.2102448,(SE006_00009:"
               "5.2329068,SE009_00009:5.2329068):0.9029422):0.2054233,"
               "SE005_00009:6.3412723):0.3714563,SE004_00009:6.7127286):"
               "0.7293362,SE003_00009:7.4420648):0.2444784,(SE002_00009:"
               "6.0534057,SE007_00009:6.0534057):1.6331375):1.594016;")
# 10 species, 9 genes (loss)
gene_tree_3 = ("(((((((SE001_02297:2.1494876,SE010_02297:2.1494876):"
               "3.7761166,SE008_02297:5.9256042):0.2102448,(SE006_02297:"
               "5.2329068,SE009_02297:5.2329068):0.9029422):0.2054233,"
               "SE005_02297:6.3412723):0.3714563,SE004_02297:6.7127286):"
               "0.7293362,SE003_02297:7.4420648):0.2444784,SE002_02297:"
               "7.6865432):1.594016;")
# MSA, 9 genes
msa_fa_3 = """>SE001/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDI\
SVGIDPVKAKKAANNRNSFSAIYKEWYEHKKQVWSVGYASELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALRTFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE002/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTQLSSIPKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWVADWLDEKLE
>SE003/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRMIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EAMQWWADWLDEKVE
>SE004/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSKITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE005/02297
MLTVKQIEKAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKFPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKNRRRCGEVFRYAIVTGAAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKGLATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMRWWADWLDEKVE
>SE006/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRL\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE008/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDNKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDNKVE
>SE009/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPMMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRL\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVNEFVFAGR\
NGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE010/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDI\
SVGIDPVKAKKASNNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
"""
# concatenated species and gene trees
species_gene_tree_1_exp = """(((((((SE001:2.1494877,SE010:1.08661):3.7761166\
,SE008:0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):0.90294\
223):0.20542323,SE005:3.0992506):0.37145632,SE004:1.8129133):0.72933621,SE00\
3:1.737411):0.24447835,(SE002:1.6606127,SE007:0.70000178):1.6331374):1.59401\
6;
(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):3.7761166,SE008_01623:5.\
9256042):0.2102448,(SE006_01623:5.2329068,SE009_01623:5.2329068):0.90294\
22):0.2054233,SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):0.7293\
362,SE003_01623:7.4420648):0.2444784,((SE002_01623:6.0534057,SE007_01623:6.0\
534057):0.4589905,((((SE001_04123:2.1494876,SE010_04123:2.1494876):3.776\
1166,SE008_04123:5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_0\
4123:5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):0.171123\
9):1.174147):1.594016;
"""
# GenBank file as input
species_genome = """
LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds, and Axl2p
            (AXL2) and Rev7p (REV7) genes, complete cds.
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;
            Saccharomycetales; Saccharomycetaceae; Saccharomyces.
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
REFERENCE   2  (bases 1 to 5028)
  AUTHORS   Roemer,T., Madden,K., Chang,J. and Snyder,M.
  TITLE     Selection of axial growth sites in yeast requires Axl2p, a novel
            plasma membrane glycoprotein
  JOURNAL   Genes Dev. 10 (7), 777-793 (1996)
  PUBMED    8846915
REFERENCE   3  (bases 1 to 5028)
  AUTHORS   Roemer,T.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-FEB-1996) Terry Roemer, Biology, Yale University, New
            Haven, CT, USA
FEATURES             Location/Qualifiers
     source          1..5028
                     /organism="Saccharomyces cerevisiae"
                     /db_xref="taxon:4932"
                     /chromosome="IX"
                     /map="9"
     CDS             <1..206
                     /codon_start=3
                     /product="TCP1-beta"
                     /protein_id="AAA98665.1"
                     /db_xref="GI:1293614"
                     /translation="SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEA
                     AEVLLRVDNIIRARPRTANRQHM"
     gene            687..3158
                     /gene="AXL2"
     CDS             687..3158
                     /gene="AXL2"
                     /note="plasma membrane glycoprotein"
                     /codon_start=1
                     /function="required for axial budding pattern of S.
                     cerevisiae"
                     /product="Axl2p"
                     /protein_id="AAA98666.1"
                     /db_xref="GI:1293615"
                     /translation="MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESF
                     TFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLYFN
                     VILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDPNE
                     VFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPE
                     TSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYV
                     YLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYG
                     DVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQ
                     DHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSA
                     NATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIA
                     CGVAIPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLN
                     NPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQ
                     SQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDS
                     YGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTK
                     HRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRL
                     VDFSNKSNVNVGQVKDIHGRIPEML"
     gene            complement(3300..4037)
                     /gene="REV7"
     CDS             complement(3300..4037)
                     /gene="REV7"
                     /codon_start=1
                     /product="Rev7p"
                     /protein_id="AAA98667.1"
                     /db_xref="GI:1293616"
                     /translation="MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQ
                     FVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVD
                     KDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDRNR
                     RVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEK
                     LISGDDKILNGVYSQYEEGESIFGSLF"
ORIGIN
        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca gtagtcagct
      121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt gcaagaccaa
      181 gaaccgccaa tagacaacat atgtaacata tttaggatat acctcgaaaa taataaaccg
      241 ccacactgtc attattataa ttagaaacag aacgcaaaaa ttatccacta tataattcaa
      301 agacgcgaaa aaaaaagaac aacgcgtcat agaacttttg gcaattcgcg tcacaaataa
      361 attttggcaa cttatgtttc ctcttcgagc agtactcgag ccctgtctca agaatgtaat
      421 aatacccatc gtaggtatgg ttaaagatag catctccaca acctcaaagc tccttgccga
      481 gagtcgccct cctttgtcga gtaattttca cttttcatat gagaacttat tttcttattc
      541 tttactctca catcctgtag tgattgacac tgcaacagcc accatcacta gaagaacaga
      601 acaattactt aatagaaaaa ttatatcttc ctcgaaacga tttcctgctt ccaacatcta
      661 cgtatatcaa gaagcattca cttaccatga cacagcttca gatttcatta ttgctgacag
      721 ctactatatc actactccat ctagtagtgg ccacgcccta tgaggcatat cctatcggaa
      781 aacaataccc cccagtggca agagtcaatg aatcgtttac atttcaaatt tccaatgata
      841 cctataaatc gtctgtagac aagacagctc aaataacata caattgcttc gacttaccga
      901 gctggctttc gtttgactct agttctagaa cgttctcagg tgaaccttct tctgacttac
      961 tatctgatgc gaacaccacg ttgtatttca atgtaatact cgagggtacg gactctgccg
     1021 acagcacgtc tttgaacaat acataccaat ttgttgttac aaaccgtcca tccatctcgc
     1081 tatcgtcaga tttcaatcta ttggcgttgt taaaaaacta tggttatact aacggcaaaa
     1141 acgctctgaa actagatcct aatgaagtct tcaacgtgac ttttgaccgt tcaatgttca
     1201 ctaacgaaga atccattgtg tcgtattacg gacgttctca gttgtataat gcgccgttac
     1261 ccaattggct gttcttcgat tctggcgagt tgaagtttac tgggacggca ccggtgataa
     1321 actcggcgat tgctccagaa acaagctaca gttttgtcat catcgctaca gacattgaag
     1381 gattttctgc cgttgaggta gaattcgaat tagtcatcgg ggctcaccag ttaactacct
     1441 ctattcaaaa tagtttgata atcaacgtta ctgacacagg taacgtttca tatgacttac
     1501 ctctaaacta tgtttatctc gatgacgatc ctatttcttc tgataaattg ggttctataa
     1561 acttattgga tgctccagac tgggtggcat tagataatgc taccatttcc gggtctgtcc
     1621 cagatgaatt actcggtaag aactccaatc ctgccaattt ttctgtgtcc atttatgata
     1681 cttatggtga tgtgatttat ttcaacttcg aagttgtctc cacaacggat ttgtttgcca
     1741 ttagttctct tcccaatatt aacgctacaa ggggtgaatg gttctcctac tattttttgc
     1801 cttctcagtt tacagactac gtgaatacaa acgtttcatt agagtttact aattcaagcc
     1861 aagaccatga ctgggtgaaa ttccaatcat ctaatttaac attagctgga gaagtgccca
     1921 agaatttcga caagctttca ttaggtttga aagcgaacca aggttcacaa tctcaagagc
     1981 tatattttaa catcattggc atggattcaa agataactca ctcaaaccac agtgcgaatg
     2041 caacgtccac aagaagttct caccactcca cctcaacaag ttcttacaca tcttctactt
     2101 acactgcaaa aatttcttct acctccgctg ctgctacttc ttctgctcca gcagcgctgc
     2161 cagcagccaa taaaacttca tctcacaata aaaaagcagt agcaattgcg tgcggtgttg
     2221 ctatcccatt aggcgttatc ctagtagctc tcatttgctt cctaatattc tggagacgca
     2281 gaagggaaaa tccagacgat gaaaacttac cgcatgctat tagtggacct gatttgaata
     2341 atcctgcaaa taaaccaaat caagaaaacg ctacaccttt gaacaacccc tttgatgatg
     2401 atgcttcctc gtacgatgat acttcaatag caagaagatt ggctgctttg aacactttga
     2461 aattggataa ccactctgcc actgaatctg atatttccag cgtggatgaa aagagagatt
     2521 ctctatcagg tatgaataca tacaatgatc agttccaatc ccaaagtaaa gaagaattat
     2581 tagcaaaacc cccagtacag cctccagaga gcccgttctt tgacccacag aataggtctt
     2641 cttctgtgta tatggatagt gaaccagcag taaataaatc ctggcgatat actggcaacc
     2701 tgtcaccagt ctctgatatt gtcagagaca gttacggatc acaaaaaact gttgatacag
     2761 aaaaactttt cgatttagaa gcaccagaga aggaaaaacg tacgtcaagg gatgtcacta
     2821 tgtcttcact ggacccttgg aacagcaata ttagcccttc tcccgtaaga aaatcagtaa
     2881 caccatcacc atataacgta acgaagcatc gtaaccgcca cttacaaaat attcaagact
     2941 ctcaaagcgg taaaaacgga atcactccca caacaatgtc aacttcatct tctgacgatt
     3001 ttgttccggt taaagatggt gaaaattttt gctgggtcca tagcatggaa ccagacagaa
     3061 gaccaagtaa gaaaaggtta gtagattttt caaataagag taatgtcaat gttggtcaag
     3121 ttaaggacat tcacggacgc atcccagaaa tgctgtgatt atacgcaacg atattttgct
     3181 taattttatt ttcctgtttt attttttatt agtggtttac agatacccta tattttattt
     3241 agtttttata cttagagaca tttaatttta attccattct tcaaatttca tttttgcact
     3301 taaaacaaag atccaaaaat gctctcgccc tcttcatatt gagaatacac tccattcaaa
     3361 attttgtcgt caccgctgat taatttttca ctaaactgat gaataatcaa aggccccacg
     3421 tcagaaccga ctaaagaagt gagttttatt ttaggaggtt gaaaaccatt attgtctggt
     3481 aaattttcat cttcttgaca tttaacccag tttgaatccc tttcaatttc tgctttttcc
     3541 tccaaactat cgaccctcct gtttctgtcc aacttatgtc ctagttccaa ttcgatcgca
     3601 ttaataactg cttcaaatgt tattgtgtca tcgttgactt taggtaattt ctccaaatgc
     3661 ataatcaaac tatttaagga agatcggaat tcgtcgaaca cttcagtttc cgtaatgatc
     3721 tgatcgtctt tatccacatg ttgtaattca ctaaaatcta aaacgtattt ttcaatgcat
     3781 aaatcgttct ttttattaat aatgcagatg gaaaatctgt aaacgtgcgt taatttagaa
     3841 agaacatcca gtataagttc ttctatatag tcaattaaag caggatgcct attaatggga
     3901 acgaactgcg gcaagttgaa tgactggtaa gtagtgtagt cgaatgactg aggtgggtat
     3961 acatttctat aaaataaaat caaattaatg tagcatttta agtataccct cagccacttc
     4021 tctacccatc tattcataaa gctgacgcaa cgattactat tttttttttc ttcttggatc
     4081 tcagtcgtcg caaaaacgta taccttcttt ttccgacctt ttttttagct ttctggaaaa
     4141 gtttatatta gttaaacagg gtctagtctt agtgtgaaag ctagtggttt cgattgactg
     4201 atattaagaa agtggaaatt aaattagtag tgtagacgta tatgcatatg tatttctcgc
     4261 ctgtttatgt ttctacgtac ttttgattta tagcaagggg aaaagaaata catactattt
     4321 tttggtaaag gtgaaagcat aatgtaaaag ctagaataaa atggacgaaa taaagagagg
     4381 cttagttcat cttttttcca aaaagcaccc aatgataata actaaaatga aaaggatttg
     4441 ccatctgtca gcaacatcag ttgtgtgagc aataataaaa tcatcacctc cgttgccttt
     4501 agcgcgtttg tcgtttgtat cttccgtaat tttagtctta tcaatgggaa tcataaattt
     4561 tccaatgaat tagcaatttc gtccaattct ttttgagctt cttcatattt gctttggaat
     4621 tcttcgcact tcttttccca ttcatctctt tcttcttcca aagcaacgat ccttctaccc
     4681 atttgctcag agttcaaatc ggcctctttc agtttatcca ttgcttcctt cagtttggct
     4741 tcactgtctt ctagctgttg ttctagatcc tggtttttct tggtgtagtt ctcattatta
     4801 gatctcaagt tattggagtc ttcagccaat tgctttgtat cagacaattg actctctaac
     4861 ttctccactt cactgtcgag ttgctcgttt ttagcggaca aagatttaat ctcgttttct
     4921 ttttcagtgt tagattgctc taattctttg agctgttctc tcagctcctc atatttttct
     4981 tgccatgact cagattctaa ttttaagcta ttcaatttct ctttgatc
//
"""


if __name__ == '__main__':
    main()
