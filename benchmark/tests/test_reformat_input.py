# -----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from os import close
from os.path import join

from skbio.util import remove_files
from skbio import TreeNode

from hgt_analysis.reformat_input import (join_trees,
                                         trim_gene_tree_leaves,
                                         species_gene_mapping)


class workflowTests(TestCase):
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

        self.files_to_remove = [self.species_tree_fp,
                                self.species_tree_2_fp,
                                self.gene_tree_1_fp,
                                self.gene_tree_2_fp,
                                self.gene_tree_3_fp,
                                self.msa_fa_3_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.working_dir)

    def test_join_trees(self):
        """ Test concatenate Newick trees into one file (species, gene)
        """
        self.output_file = join(self.working_dir, 'output_file.nwk')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        join_trees(gene_tree_1, species_tree, self.output_file)
        with open(self.output_file, 'U') as out_f:
            species_gene_tree_1_obs = out_f.read()
        self.assertEqual(species_gene_tree_1_obs, species_gene_tree_1_exp)

    def test_trim_gene_tree_leaves(self):
        """ Test remove '_GENENAME' from tree leaf names (if exists)
        """
        leaves_exp = ["SE001", "SE002", "SE003", "SE004", "SE005", "SE006",
                      "SE007", "SE008", "SE009", "SE010"]
        gene_tree_2 = TreeNode.read(self.gene_tree_2_fp, format='newick')
        trim_gene_tree_leaves(gene_tree_2)
        leaves_obs = []
        for node in gene_tree_2.tips():
            leaves_obs.append(node.name)
        self.assertItemsEqual(leaves_obs, leaves_exp)

    def test_species_gene_mapping(self):
        """ Test finding the association between species and gene tree leaves
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        mapping_exp = {"SE001": ["SE001_01623", "SE001_04123"],
                       "SE002": ["SE002_01623"], "SE003": ["SE003_01623"],
                       "SE004": ["SE004_01623"],
                       "SE005": ["SE005_04123", "SE005_01623"],
                       "SE006": ["SE006_01623", "SE006_04123"],
                       "SE007": ["SE007_01623"],
                       "SE008": ["SE008_01623", "SE008_04123"],
                       "SE009": ["SE009_01623", "SE009_04123"],
                       "SE010": ["SE010_04123", "SE010_01623"]}
        mapping_obs = species_gene_mapping(gene_tree_1, species_tree)
        self.assertItemsEqual(mapping_obs, mapping_exp)

    def test_species_gene_mapping_check_species_labels(self):
        species_tree = TreeNode.read(self.species_tree_2_fp, format='newick')
        gene_tree_3 = TreeNode.read(self.gene_tree_3_fp, format='newick')
        self.assertRaises(ValueError,
                          species_gene_mapping,
                          gene_tree=gene_tree_3,
                          species_tree=species_tree)


# 10 species
species_tree = """(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:1.6606127,SE007:0.70000178):1.6331374):1.594016;"""
# 10 species, 1 with same label
species_tree_2 = """(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:0.86305436):0.21024487,(SE001:0.56704221,SE009:0.5014676):0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:1.6606127,SE007:0.70000178):1.6331374):1.594016;"""
# 10 species, 16 genes (gain)
gene_tree_1 = """(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):3.7761166,SE008_01623:5.9256042):0.2102448,(SE006_01623:5.2329068,SE009_01623:5.2329068):0.9029422):0.2054233,SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):0.7293362,SE003_01623:7.4420648):0.2444784,((SE002_01623:6.0534057,SE007_01623:6.0534057):0.4589905,((((SE001_04123:2.1494876,SE010_04123:2.1494876):3.7761166,SE008_04123:5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_04123:5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):0.1711239):1.174147):1.594016;"""
# 10 species, 10 genes
gene_tree_2 = """(((((((SE001_00009:2.1494876,SE010_00009:2.1494876):3.7761166,SE008_00009:5.9256042):0.2102448,(SE006_00009:5.2329068,SE009_00009:5.2329068):0.9029422):0.2054233,SE005_00009:6.3412723):0.3714563,SE004_00009:6.7127286):0.7293362,SE003_00009:7.4420648):0.2444784,(SE002_00009:6.0534057,SE007_00009:6.0534057):1.6331375):1.594016;"""
# 10 species, 9 genes (loss)
gene_tree_3 = """(((((((SE001_02297:2.1494876,SE010_02297:2.1494876):3.7761166,SE008_02297:5.9256042):0.2102448,(SE006_02297:5.2329068,SE009_02297:5.2329068):0.9029422):0.2054233,SE005_02297:6.3412723):0.3714563,SE004_02297:6.7127286):0.7293362,SE003_02297:7.4420648):0.2444784,SE002_02297:7.6865432):1.594016;"""
# MSA, 9 genes
msa_fa_3 = """>SE001/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDISVGIDPVKAKKAANNRNSFSAIYKEWYEHKKQVWSVGYASELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALRTFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDEKVE
>SE002/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTQLSSIPKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWVADWLDEKLE
>SE003/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRMIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRREAMQWWADWLDEKVE
>SE004/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSKITKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDEKVE
>SE005/02297
MLTVKQIEKAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKFPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMERANKNRRRCGEVFRYAIVTGAAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKGLATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMRWWADWLDEKVE
>SE006/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRLEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDEKVE
>SE008/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNDNKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDNKVE
>SE009/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPMMTLQEARDKAWTARKDISVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRLEDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVNEFVFAGRNGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDEKVE
>SE010/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDISVGIDPVKAKKASNNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRFEDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALATFSGSIVSLIATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGRNDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRREMMQWWADWLDEKVE
"""
# concatenated species and gene trees
species_gene_tree_1_exp = """(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:1.6606127,SE007:0.70000178):1.6331374):1.594016;
(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):3.7761166,SE008_01623:5.9256042):0.2102448,(SE006_01623:5.2329068,SE009_01623:5.2329068):0.9029422):0.2054233,SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):0.7293362,SE003_01623:7.4420648):0.2444784,((SE002_01623:6.0534057,SE007_01623:6.0534057):0.4589905,((((SE001_04123:2.1494876,SE010_04123:2.1494876):3.7761166,SE008_04123:5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_04123:5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):0.1711239):1.174147):1.594016;
"""


if __name__ == '__main__':
    main()
