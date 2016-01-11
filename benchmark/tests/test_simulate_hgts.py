# ----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# dependencies: scikit-bio development

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os import makedirs
from os.path import join, exists, dirname, abspath
import time
import copy
from operator import itemgetter

from skbio import Sequence, SequenceCollection

from simulate_hgts import (extract_genbank,
                           launch_orthofinder,
                           parse_orthofinder,
                           simulate_orthologous_rep,
                           simulate_novel_acq,
                           write_results,
                           simulate_genbank,
                           simulate_azad_lawrence)


class SimulateHGTsTests(TestCase):
    """Tests for simulating HGTs in artificial and genuine genomes
    """

    def setUp(self):
        """Set up working directory and test files
        """
        # Candidatus Carsonella ruddii PV and Candidatus Carsonella ruddii DC
        # genomes (~200 genes each)
        self.root = join(dirname(abspath(__file__)), "data")
        self.working_dir = mkdtemp()
        self.simulated_dir = join(self.working_dir, "simulated")
        makedirs(self.simulated_dir)
        self.proteomes_dir = join(self.working_dir, "proteomes")
        makedirs(self.proteomes_dir)

        # seqs 1
        self.seqs_prot_1_fp = join(self.proteomes_dir, "seqs_prot_1.fasta")
        with open(self.seqs_prot_1_fp, 'w') as tmp:
            tmp.write(seqs_prot_1)
        # seqs 2
        self.seqs_prot_2_fp = join(self.proteomes_dir, "seqs_prot_2.fasta")
        with open(self.seqs_prot_2_fp, 'w') as tmp:
            tmp.write(seqs_prot_2)

        self.genes_donor = {
                       'D_1': ['MEKEYSPKKIENYVQEFWKK', 30, 90, '+'],
                       'D_2': ['MKKNIILNLIGLRCPEPIMI', 140, 200, '+'],
                       'D_3': ['MNKILVIILFSLVSLTWGTTWIAMKIA', 220, 301, '-']}
        self.seq_donor = Sequence(
                      "TTACAGAATTTGTAGATCCATGTATTTTTGATGGAGAAGGAGTACAGCCCCAA"
                      "GAAGATCGAGAACTACGTGCAGGAGTTCTGGAAGAAGTTTATATATCTTTTTT"
                      "TATTAATAATTAATATAAAAAAATTTTCTATATTATGAAGAAGAACATCATCC"
                      "TGAACCTGATCGGCCTGAGGTGCCCCGAGCCCATCATGATCCGTAGTCGTGAT"
                      "GCTGATGTATGAACAAGATCCTGGTGATCATCCTGTTCAGCCTGGTGAGCCTG"
                      "ACCTGGGGCACCACCTGGATCGCCATGAAGATCGCCTGTATATACCAGCAATT"
                      "TCCCCAATTTTGCTTTTAAAATTTGAAATTGATTTTTTTATTTTAGAAAACGT"
                      "TGGTTTTTGACCAGTAATATATTTTATTGA")
        self.genes_recip = {
                       'R_1': ['MNLEYNPKKIESFVQQYWRN', 45, 104, '+'],
                       'R_2': ['MTELITLNLLGLRCPEPLMV', 120, 179, '+'],
                       'R_3': ['MWGTTWIAMKIVITTIPPIFATGLRFL', 240, 320, '+']}
        self.seq_recip = Sequence(
                      "AATTAGCCTATTAAATTTATTAATTTTTATATTGCTTAATATATAATGAATTT"
                      "GGAATATAATCCAAAAAAAATTGAATCTTTTGTTCAACAATATTGGAGAAATT"
                      "TAATATTAATTTGAATGACTGAATTGATTACTTTGAATTTGTTGGGTTTGAGA"
                      "TGTCCAGAACCATTGATGGTTTTTTTTTACTAATTTTTATTATTAACAAAAAA"
                      "TTTAAAGATATTTTTAAAAAATTCAATAATGTGGGGTACTACTTGGATTGCTA"
                      "TGAAAATTGTTATTACTACTATTCCACCAATTTTTGCTACTGGTTTGAGATTT"
                      "TTGATTTTAGACTGGTATTTCAAGAACGATTACTTTTAAACTGGCGTTTAAAT"
                      "ATCAACATCTCCCAGCTATCCTACACAAAAA")

    def tearDown(self):
        rmtree(self.working_dir)

    def test_extract_genbank(self):
        """Test parsing sequence and gene information from a GenBank record.
        """
        seq, genes = extract_genbank(
            join(self.root, "genbank_sample_record.gbk"))
        genes_exp = {'AAA98665.1': ['SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKR'
                                    'AVVSSASEAAEVLLRVDNIIRARPRTANRQHM', 0,
                                    205, '+'],
                     'AAA98667.1': ['MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYT'
                                    'TYQSFNLPQFVPINRHPALIDYIEELILDVLSKLT'
                                    'HVYRFSICIINKKNDLCIEKYVLDFSELQHVDKDD'
                                    'QIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTIT'
                                    'FEAVINAIELELGHKLDRNRRVDSLEEKAEIERDS'
                                    'NWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPL'
                                    'IIHQFSEKLISGDDKILNGVYSQYEEGESIFGSLF',
                                    3299, 4036, '-'],
                     'AAA98666.1': ['MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYP'
                                    'PVARVNESFTFQISNDTYKSSVDKTAQITYNCFDL'
                                    'PSWLSFDSSSRTFSGEPSSDLLSDANTTLYFNVIL'
                                    'EGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLA'
                                    'LLKNYGYTNGKNALKLDPNEVFNVTFDRSMFTNEE'
                                    'SIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPV'
                                    'INSAIAPETSYSFVIIATDIEGFSAVEVEFELVIG'
                                    'AHQLTTSIQNSLIINVTDTGNVSYDLPLNYVYLDD'
                                    'DPISSDKLGSINLLDAPDWVALDNATISGSVPDEL'
                                    'LGKNSNPANFSVSIYDTYGDVIYFNFEVVSTTDLF'
                                    'AISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSL'
                                    'EFTNSSQDHDWVKFQSSNLTLAGEVPKNFDKLSLG'
                                    'LKANQGSQSQELYFNIIGMDSKITHSNHSANATST'
                                    'RSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAA'
                                    'LPAANKTSSHNKKAVAIACGVAIPLGVILVALICF'
                                    'LIFWRRRRENPDDENLPHAISGPDLNNPANKPNQE'
                                    'NATPLNNPFDDDASSYDDTSIARRLAALNTLKLDN'
                                    'HSATESDISSVDEKRDSLSGMNTYNDQFQSQSKEE'
                                    'LLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKS'
                                    'WRYTGNLSPVSDIVRDSYGSQKTVDTEKLFDLEAP'
                                    'EKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSP'
                                    'YNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSSSD'
                                    'DFVPVKDGENFCWVHSMEPDRRPSKKRLVDFSNKS'
                                    'NVNVGQVKDIHGRIPEML',
                                    686, 3157, '+']}
        self.assertEqual(str(seq), sample_seq)
        self.assertDictEqual(genes, genes_exp)

    def test_launch_orthofinder(self):
        """Test running OrthoFinder.
        """
        launch_orthofinder(self.proteomes_dir, 1)
        date = time.strftime("%c").split()
        results_dir = join(
            self.proteomes_dir, "Results_%s%s" % (date[1], date[2]))
        orthogroups_exp = [['YP_002468181.1', 'YP_004590122.1'],
                           ['YP_004590123.1', 'YP_002468184.1'],
                           ['YP_002468032.1', 'YP_004590028.1']]
        orthogroups_act = []
        with open(join(results_dir, "OrthologousGroups.txt"), 'U') as o:
            orthogroups_act = [line.split()[1:] for line in o]

        self.assertItemsEqual(orthogroups_act, orthogroups_exp)

    def test_parse_orthofinder(self):
        """Test parsing OrthoFinder results.
        """
        launch_orthofinder(self.proteomes_dir, 1)
        date = time.strftime("%c").split()
        results_dir = join(
            self.proteomes_dir, "Results_%s%s" % (date[1], date[2]),
            "WorkingDirectory")
        species_ids, sequence_ids, orthologous_groups =\
            parse_orthofinder(results_dir)
        species_ids_exp = {'1': 'seqs_prot_2.fasta', '0': 'seqs_prot_1.fasta'}
        seq_ids_exp = {'1_2': 'YP_004590028.1', '1_1': 'YP_004590123.1',
                       '1_0': 'YP_004590122.1', '0_2': 'YP_002468032.1',
                       '0_0': 'YP_002468181.1', '0_1': 'YP_002468184.1'}
        orthogroups_exp = [['0_0', '1_0'], ['0_1', '1_1'], ['0_2', '1_2']]
        self.assertDictEqual(species_ids, species_ids_exp)
        self.assertDictEqual(sequence_ids, seq_ids_exp)
        self.assertItemsEqual(orthologous_groups, orthogroups_exp)

    def test_simulate_orthologous_rep(self):
        """Test simulating orthologous gene replacement HGTs.
        """
        # genes_recip and seq_recip are modified within function
        # simulate_orthologous_rep, store their original copies for
        # downstream verification
        genes_donor_orig = copy.deepcopy(self.genes_donor)
        genes_recip_orig = copy.deepcopy(self.genes_recip)
        sequence_ids = {'0_0': 'D_1', '0_1': 'D_2', '0_2': 'D_3',
                        '1_0': 'R_1', '1_1': 'R_2', '1_2': 'R_3'}
        orthologous_groups = [['0_0', '1_0'], ['0_1', '1_1'], ['0_2', '1_2']]
        orthologous_rep_prob = 0.5
        percentage_hgts = 0.5
        log_fp = join(self.working_dir, "log.txt")
        with open(log_fp, 'w') as log_f:
            seq_recip = simulate_orthologous_rep(self.genes_donor,
                                                 self.seq_donor,
                                                 self.genes_recip,
                                                 self.seq_recip,
                                                 sequence_ids,
                                                 orthologous_groups,
                                                 orthologous_rep_prob,
                                                 percentage_hgts,
                                                 log_f)
        # number of genes in donor and recipient genome should remain the
        # same
        self.assertEqual(len(self.genes_recip), len(genes_recip_orig))
        self.assertItemsEqual(self.genes_donor, genes_donor_orig)
        translated_nucl = {'D_1': "ATGGAGAAGGAGTACAGCCCCAAGAAGATCGAGAACTACGTG"
                                  "CAGGAGTTCTGGAAGAAG",
                           'D_2': "ATGAAGAAGAACATCATCCTGAACCTGATCGGCCTGAGGTGC"
                                  "CCCGAGCCCATCATGATC",
                           'D_3': "ATGAACAAGATCCTGGTGATCATCCTGTTCAGCCTGGTGAGC"
                                  "CTGACCTGGGGCACCACCTGGATCGCCATGAAGATCGCC",
                           'R_1': "ATGAATTTGGAATATAATCCAAAAAAAATTGAATCTTTTGTT"
                                  "CAACAATATTGGAGAAA",
                           'R_2': "ATGACTGAATTGATTACTTTGAATTTGTTGGGTTTGAGATGT"
                                  "CCAGAACCATTGATGGT",
                           'R_3': "ATGTGGGGTACTACTTGGATTGCTATGAAAATTGTTATTACT"
                                  "ACTATTCCACCAATTTTTGCTACTGGTTTGAGATTTTT"}
        with open(log_fp, 'U') as log_f:
            hgts_sim = [line.strip().split()
                        for line in log_f
                        if not line.startswith('#')]
        # verify HGTs
        for hgt in hgts_sim:
            recip_label_replaced = hgt[4]
            self.assertTrue(recip_label_replaced not in self.genes_recip)
            hgt_gene = hgt[5]
            self.assertTrue(hgt_gene in self.genes_recip)
            donor_label = hgt[1]
            # protein sequences are equal for donor and recipient gene
            self.assertEqual(genes_donor_orig[donor_label][0],
                             self.genes_recip[hgt_gene][0])
            recip_start = hgt[6]
            recip_end = hgt[7]
            # length of gene (nucleotide format) is consistent with the
            # gene start and end positions on the recipient genome
            self.assertEqual(len(self.genes_recip[hgt_gene][0])*3,
                             int(recip_end)-int(recip_start))
            # genome subsequence representing HGT in recipient genome is
            # equal to the nucleotide format of the donor gene
            self.assertEqual(
                str(seq_recip[
                    self.genes_recip[hgt_gene][1]:self.genes_recip[hgt_gene][2]]),
                    translated_nucl[donor_label])

    def test_simulate_novel_acq(self):
        """Test simulating novel gene acquisition HGTs.
        """
               # genes_recip and seq_recip are modified within function
        # simulate_orthologous_rep, store their original copies for
        # downstream verification
        genes_donor_orig = copy.deepcopy(self.genes_donor)
        genes_recip_orig = copy.deepcopy(self.genes_recip)
        sequence_ids = {'0_0': 'D_1', '0_1': 'D_2', '0_2': 'D_3',
                        '1_0': 'R_1', '1_1': 'R_2', '1_2': 'R_3'}
        orthologous_groups = [['0_0', '1_0'], ['0_1', '1_1'], ['0_2', '1_2']]
        orthologous_rep_prob = 0.0
        percentage_hgts = 0.5
        log_fp = join(self.working_dir, "log.txt")
        with open(log_fp, 'w') as log_f:
            seq_recip = simulate_novel_acq(self.genes_donor,
                                           self.seq_donor,
                                           self.genes_recip,
                                           self.seq_recip,
                                           orthologous_rep_prob,
                                           percentage_hgts,
                                           log_f)
        # number of genes in donor genome should remain the same
        self.assertItemsEqual(self.genes_donor, genes_donor_orig)
        translated_nucl = {'D_1': "ATGGAGAAGGAGTACAGCCCCAAGAAGATCGAGAACTACGTG"
                                  "CAGGAGTTCTGGAAGAAG",
                           'D_2': "ATGAAGAAGAACATCATCCTGAACCTGATCGGCCTGAGGTGC"
                                  "CCCGAGCCCATCATGATC",
                           'D_3': "ATGAACAAGATCCTGGTGATCATCCTGTTCAGCCTGGTGAGC"
                                  "CTGACCTGGGGCACCACCTGGATCGCCATGAAGATCGCC",
                           'R_1': "ATGAATTTGGAATATAATCCAAAAAAAATTGAATCTTTTGTT"
                                  "CAACAATATTGGAGAAA",
                           'R_2': "ATGACTGAATTGATTACTTTGAATTTGTTGGGTTTGAGATGT"
                                  "CCAGAACCATTGATGGT",
                           'R_3': "ATGTGGGGTACTACTTGGATTGCTATGAAAATTGTTATTACT"
                                  "ACTATTCCACCAATTTTTGCTACTGGTTTGAGATTTTT"}
        with open(log_fp, 'U') as log_f:
            hgts_sim = [line.strip().split()
                        for line in log_f
                        if not line.startswith('#')]
        gene_positions =\
            [(self.genes_recip[g][1], self.genes_recip[g][2], True)
             if 'hgt' in g
             else (self.genes_recip[g][1], self.genes_recip[g][2], False)
             for g in self.genes_recip]
        # verify HGTs
        hgt_genes = {}
        for hgt in hgts_sim:
            donor_label = hgt[1]
            hgt_label = hgt[4]
            recip_start = int(hgt[5])
            recip_end = int(hgt[6])
            self.assertTrue(hgt_label in self.genes_recip)
            self.assertEqual(recip_start, self.genes_recip[hgt_label][1])
            self.assertEqual(recip_end, self.genes_recip[hgt_label][2])
            self.assertEqual(hgt[7], self.genes_recip[hgt_label][3])
            # length of gene (nucleotide format) is consistent with the
            # gene start and end positions on the recipient genome
            self.assertEqual(len(self.genes_recip[hgt_label][0])*3,
                             int(recip_end)-int(recip_start))
            # genome subsequence representing HGT in recipient genome is
            # equal to the nucleotide format of the donor gene
            self.assertEqual(
                str(seq_recip[
                    self.genes_recip[hgt_label][1]:self.genes_recip[hgt_label][2]]),
                    translated_nucl[donor_label])
        gene_positions_s = sorted(gene_positions, key=itemgetter(0))
        # verify none of the original genes overlap with HGTs
        for x in xrange(0, len(gene_positions_s)):
            # HGT gene
            if gene_positions_s[x][2]:
                if x < len(gene_positions_s)-1:
                    self.assertTrue(
                        gene_positions_s[x][1] < gene_positions_s[x+1][0])
                if x > 0:
                    self.assertTrue(
                        gene_positions_s[x][0] > gene_positions_s[x-1][1])

    def test_write_results(self):
        """Test writing HGT results to FASTA files.
        """
        self.genes_recip['D_2_hgt_n'] = ['MKKNIILNLIGLRCPEPIMI', 321, 381, '+']
        donor_genbank_fp = join(self.proteomes_dir, "donor.fna")
        recipient_genbank_fp = join(self.proteomes_dir, "recip.fna")
        dnr_g_nucl_fp, dnr_g_aa_fp, rcp_g_nucl_fp, rcp_g_aa_fp =\
            write_results(self.genes_donor,
                          donor_genbank_fp,
                          self.genes_recip,
                          recipient_genbank_fp,
                          self.seq_donor,
                          self.seq_recip,
                          self.simulated_dir)
        donor_nucl = Sequence.read(dnr_g_nucl_fp, format='fasta')
        # test for correctness of donor nucleotide genome sequence
        self.assertEqual(str(donor_nucl), str(self.seq_donor))
        recip_nucl = Sequence.read(rcp_g_nucl_fp, format='fasta')
        # test for correctness of recipient nucleotide genome sequence
        self.assertEqual(str(recip_nucl), str(self.seq_recip))
        donor_aa = SequenceCollection.read(dnr_g_aa_fp, format='fasta')
        donor_aa_dict = {}
        for seq in donor_aa:
            donor_aa_dict[seq.metadata['id']] = seq
        # test for correctness of donor protein coding sequences
        self.assertItemsEqual(donor_aa_dict, self.genes_donor)
        recip_aa = SequenceCollection.read(rcp_g_aa_fp, format='fasta')
        recip_aa_dict = {}
        for seq in recip_aa:
            recip_aa_dict[seq.metadata['id']] = seq
        # test for correctness of recipient protein coding sequences
        self.assertItemsEqual(recip_aa_dict, self.genes_recip)

    def test_simulate_genbank(self):
        """Test simulating HGTs using input donor and recipient GenBank files.
        """
        donor_genbank_fp = join(self.root, "genbank",
                                "GCF_000010365.1_ASM1036v1_genomic.gbff")
        recip_genbank_fp = join(self.root, "genbank",
                                "GCF_000441575.1_ASM44157v1_genomic.gbff")
        output_dir = self.simulated_dir
        percentage_hgts = 0.2
        orthologous_rep_prob = 0.5
        log_fp = join(self.working_dir, "log.txt")
        threads = 1
        with open(log_fp, 'w') as log_f:
            dnr_nucl_fp, dnr_aa_fp, rcp_nucl_fp, rcp_aa_fp =\
                simulate_genbank(donor_genbank_fp,
                                 recip_genbank_fp,
                                 output_dir,
                                 percentage_hgts,
                                 orthologous_rep_prob,
                                 log_f,
                                 threads)
        with open(log_fp, 'U') as log_f:
            hgts_sim = [line.strip().split()[4] if line[0]=='n' else
                        line.strip().split()[5]
                        for line in log_f
                        if not line.startswith('#')]
        donor_aa = SequenceCollection.read(dnr_aa_fp, format='fasta')
        rcp_aa = SequenceCollection.read(rcp_aa_fp, format='fasta')
        rcp_hgts = [_id for _id in rcp_aa.ids() if 'hgt' in _id]
        self.assertItemsEqual(rcp_hgts, hgts_sim)

    def test_simulate_azad_lawrence(test):
        """Test simulating HGTs using artificial genomes by Azad et al., 2005.
        """
        donor_artificial_fp = join(self.root, "genbank",
                                "GCF_000010365.1_ASM1036v1_genomic.gbff")
        donor_artificial_annotation_fp = join(self.root, "genbank",
                                "GCF_000441575.1_ASM44157v1_genomic.gbff")

simulate_azad_lawrence(donor_artificial_fp,
                           donor_artificial_annotation_fp,
                           recip_artificial_fp,
                           recip_artificial_annotation_fp,
                           output_dir,
                           percentage_hgts,
                           orthologous_rep_prob,
                           log_f,
                           threads,
                           verbose=False)

sample_seq = """GATCCTCCATATACAACGGTATCTCCACCTCAGGTTTAGATCTCAACAACGGAACCATTGC\
CGACATGAGACAGTTAGGTATCGTCGAGAGTTACAAGCTAAAACGAGCAGTAGTCAGCTCTGCATCTGAAGCCGCTG\
AAGTTCTACTAAGGGTGGATAACATCATCCGTGCAAGACCAAGAACCGCCAATAGACAACATATGTAACATATTTAG\
GATATACCTCGAAAATAATAAACCGCCACACTGTCATTATTATAATTAGAAACAGAACGCAAAAATTATCCACTATA\
TAATTCAAAGACGCGAAAAAAAAAGAACAACGCGTCATAGAACTTTTGGCAATTCGCGTCACAAATAAATTTTGGCA\
ACTTATGTTTCCTCTTCGAGCAGTACTCGAGCCCTGTCTCAAGAATGTAATAATACCCATCGTAGGTATGGTTAAAG\
ATAGCATCTCCACAACCTCAAAGCTCCTTGCCGAGAGTCGCCCTCCTTTGTCGAGTAATTTTCACTTTTCATATGAG\
AACTTATTTTCTTATTCTTTACTCTCACATCCTGTAGTGATTGACACTGCAACAGCCACCATCACTAGAAGAACAGA\
ACAATTACTTAATAGAAAAATTATATCTTCCTCGAAACGATTTCCTGCTTCCAACATCTACGTATATCAAGAAGCAT\
TCACTTACCATGACACAGCTTCAGATTTCATTATTGCTGACAGCTACTATATCACTACTCCATCTAGTAGTGGCCAC\
GCCCTATGAGGCATATCCTATCGGAAAACAATACCCCCCAGTGGCAAGAGTCAATGAATCGTTTACATTTCAAATTT\
CCAATGATACCTATAAATCGTCTGTAGACAAGACAGCTCAAATAACATACAATTGCTTCGACTTACCGAGCTGGCTT\
TCGTTTGACTCTAGTTCTAGAACGTTCTCAGGTGAACCTTCTTCTGACTTACTATCTGATGCGAACACCACGTTGTA\
TTTCAATGTAATACTCGAGGGTACGGACTCTGCCGACAGCACGTCTTTGAACAATACATACCAATTTGTTGTTACAA\
ACCGTCCATCCATCTCGCTATCGTCAGATTTCAATCTATTGGCGTTGTTAAAAAACTATGGTTATACTAACGGCAAA\
AACGCTCTGAAACTAGATCCTAATGAAGTCTTCAACGTGACTTTTGACCGTTCAATGTTCACTAACGAAGAATCCAT\
TGTGTCGTATTACGGACGTTCTCAGTTGTATAATGCGCCGTTACCCAATTGGCTGTTCTTCGATTCTGGCGAGTTGA\
AGTTTACTGGGACGGCACCGGTGATAAACTCGGCGATTGCTCCAGAAACAAGCTACAGTTTTGTCATCATCGCTACA\
GACATTGAAGGATTTTCTGCCGTTGAGGTAGAATTCGAATTAGTCATCGGGGCTCACCAGTTAACTACCTCTATTCA\
AAATAGTTTGATAATCAACGTTACTGACACAGGTAACGTTTCATATGACTTACCTCTAAACTATGTTTATCTCGATG\
ACGATCCTATTTCTTCTGATAAATTGGGTTCTATAAACTTATTGGATGCTCCAGACTGGGTGGCATTAGATAATGCT\
ACCATTTCCGGGTCTGTCCCAGATGAATTACTCGGTAAGAACTCCAATCCTGCCAATTTTTCTGTGTCCATTTATGA\
TACTTATGGTGATGTGATTTATTTCAACTTCGAAGTTGTCTCCACAACGGATTTGTTTGCCATTAGTTCTCTTCCCA\
ATATTAACGCTACAAGGGGTGAATGGTTCTCCTACTATTTTTTGCCTTCTCAGTTTACAGACTACGTGAATACAAAC\
GTTTCATTAGAGTTTACTAATTCAAGCCAAGACCATGACTGGGTGAAATTCCAATCATCTAATTTAACATTAGCTGG\
AGAAGTGCCCAAGAATTTCGACAAGCTTTCATTAGGTTTGAAAGCGAACCAAGGTTCACAATCTCAAGAGCTATATT\
TTAACATCATTGGCATGGATTCAAAGATAACTCACTCAAACCACAGTGCGAATGCAACGTCCACAAGAAGTTCTCAC\
CACTCCACCTCAACAAGTTCTTACACATCTTCTACTTACACTGCAAAAATTTCTTCTACCTCCGCTGCTGCTACTTC\
TTCTGCTCCAGCAGCGCTGCCAGCAGCCAATAAAACTTCATCTCACAATAAAAAAGCAGTAGCAATTGCGTGCGGTG\
TTGCTATCCCATTAGGCGTTATCCTAGTAGCTCTCATTTGCTTCCTAATATTCTGGAGACGCAGAAGGGAAAATCCA\
GACGATGAAAACTTACCGCATGCTATTAGTGGACCTGATTTGAATAATCCTGCAAATAAACCAAATCAAGAAAACGC\
TACACCTTTGAACAACCCCTTTGATGATGATGCTTCCTCGTACGATGATACTTCAATAGCAAGAAGATTGGCTGCTT\
TGAACACTTTGAAATTGGATAACCACTCTGCCACTGAATCTGATATTTCCAGCGTGGATGAAAAGAGAGATTCTCTA\
TCAGGTATGAATACATACAATGATCAGTTCCAATCCCAAAGTAAAGAAGAATTATTAGCAAAACCCCCAGTACAGCC\
TCCAGAGAGCCCGTTCTTTGACCCACAGAATAGGTCTTCTTCTGTGTATATGGATAGTGAACCAGCAGTAAATAAAT\
CCTGGCGATATACTGGCAACCTGTCACCAGTCTCTGATATTGTCAGAGACAGTTACGGATCACAAAAAACTGTTGAT\
ACAGAAAAACTTTTCGATTTAGAAGCACCAGAGAAGGAAAAACGTACGTCAAGGGATGTCACTATGTCTTCACTGGA\
CCCTTGGAACAGCAATATTAGCCCTTCTCCCGTAAGAAAATCAGTAACACCATCACCATATAACGTAACGAAGCATC\
GTAACCGCCACTTACAAAATATTCAAGACTCTCAAAGCGGTAAAAACGGAATCACTCCCACAACAATGTCAACTTCA\
TCTTCTGACGATTTTGTTCCGGTTAAAGATGGTGAAAATTTTTGCTGGGTCCATAGCATGGAACCAGACAGAAGACC\
AAGTAAGAAAAGGTTAGTAGATTTTTCAAATAAGAGTAATGTCAATGTTGGTCAAGTTAAGGACATTCACGGACGCA\
TCCCAGAAATGCTGTGATTATACGCAACGATATTTTGCTTAATTTTATTTTCCTGTTTTATTTTTTATTAGTGGTTT\
ACAGATACCCTATATTTTATTTAGTTTTTATACTTAGAGACATTTAATTTTAATTCCATTCTTCAAATTTCATTTTT\
GCACTTAAAACAAAGATCCAAAAATGCTCTCGCCCTCTTCATATTGAGAATACACTCCATTCAAAATTTTGTCGTCA\
CCGCTGATTAATTTTTCACTAAACTGATGAATAATCAAAGGCCCCACGTCAGAACCGACTAAAGAAGTGAGTTTTAT\
TTTAGGAGGTTGAAAACCATTATTGTCTGGTAAATTTTCATCTTCTTGACATTTAACCCAGTTTGAATCCCTTTCAA\
TTTCTGCTTTTTCCTCCAAACTATCGACCCTCCTGTTTCTGTCCAACTTATGTCCTAGTTCCAATTCGATCGCATTA\
ATAACTGCTTCAAATGTTATTGTGTCATCGTTGACTTTAGGTAATTTCTCCAAATGCATAATCAAACTATTTAAGGA\
AGATCGGAATTCGTCGAACACTTCAGTTTCCGTAATGATCTGATCGTCTTTATCCACATGTTGTAATTCACTAAAAT\
CTAAAACGTATTTTTCAATGCATAAATCGTTCTTTTTATTAATAATGCAGATGGAAAATCTGTAAACGTGCGTTAAT\
TTAGAAAGAACATCCAGTATAAGTTCTTCTATATAGTCAATTAAAGCAGGATGCCTATTAATGGGAACGAACTGCGG\
CAAGTTGAATGACTGGTAAGTAGTGTAGTCGAATGACTGAGGTGGGTATACATTTCTATAAAATAAAATCAAATTAA\
TGTAGCATTTTAAGTATACCCTCAGCCACTTCTCTACCCATCTATTCATAAAGCTGACGCAACGATTACTATTTTTT\
TTTTCTTCTTGGATCTCAGTCGTCGCAAAAACGTATACCTTCTTTTTCCGACCTTTTTTTTAGCTTTCTGGAAAAGT\
TTATATTAGTTAAACAGGGTCTAGTCTTAGTGTGAAAGCTAGTGGTTTCGATTGACTGATATTAAGAAAGTGGAAAT\
TAAATTAGTAGTGTAGACGTATATGCATATGTATTTCTCGCCTGTTTATGTTTCTACGTACTTTTGATTTATAGCAA\
GGGGAAAAGAAATACATACTATTTTTTGGTAAAGGTGAAAGCATAATGTAAAAGCTAGAATAAAATGGACGAAATAA\
AGAGAGGCTTAGTTCATCTTTTTTCCAAAAAGCACCCAATGATAATAACTAAAATGAAAAGGATTTGCCATCTGTCA\
GCAACATCAGTTGTGTGAGCAATAATAAAATCATCACCTCCGTTGCCTTTAGCGCGTTTGTCGTTTGTATCTTCCGT\
AATTTTAGTCTTATCAATGGGAATCATAAATTTTCCAATGAATTAGCAATTTCGTCCAATTCTTTTTGAGCTTCTTC\
ATATTTGCTTTGGAATTCTTCGCACTTCTTTTCCCATTCATCTCTTTCTTCTTCCAAAGCAACGATCCTTCTACCCA\
TTTGCTCAGAGTTCAAATCGGCCTCTTTCAGTTTATCCATTGCTTCCTTCAGTTTGGCTTCACTGTCTTCTAGCTGT\
TGTTCTAGATCCTGGTTTTTCTTGGTGTAGTTCTCATTATTAGATCTCAAGTTATTGGAGTCTTCAGCCAATTGCTT\
TGTATCAGACAATTGACTCTCTAACTTCTCCACTTCACTGTCGAGTTGCTCGTTTTTAGCGGACAAAGATTTAATCT\
CGTTTTCTTTTTCAGTGTTAGATTGCTCTAATTCTTTGAGCTGTTCTCTCAGCTCCTCATATTTTTCTTGCCATGAC\
TCAGATTCTAATTTTAAGCTATTCAATTTCTCTTTGATC\
"""

seqs_prot_1 = """>YP_002468181.1
MEKEYSPKKIENYVQEFWKKNKTFEVKEDPKKEKYYCLPMLPYPSGKLHMGHVRNYTISDVISRYQRMLGKNVLQPM\
GWDAFGLPAEEAAIRNNTDPFSWTQKNIKYMKKQLQSLGFSYDWSREITTCHPEYYHWEQWFFTKLYEKKLVYKKNS\
LVNWCSYDKTVLANEQVIDGCCWRCQNKIRIKQIPQWFIKIRNYAESLYQDLKKLTHWPENVKNMQRNWIGRIKGFE\
ITLNVFNTCQKLKVFTQRLDLLMGVTYISISSCHKLSINLSKKNELIKKFIKKYRYISQEEQYKVKYTGINTNLFVV\
HPITKKTIPIWISNATHIEYGTNAVLSIPGHNENDWNFAVKNNLKIKYVIFNPDHQEPNLYTSFLDIKGTLFNSQEF\
NGLNLKDGTEKIKKILYKKKILKEKINYKLQDWCISRQRYWGTPIPMAKFKNGKMIPIPENQLPVVLPKIRKNTNLL\
QQAINFNSKWAEIFIHGKHAIREIDTFDTFMESSWYYARYTCPNFNTGMIDSIASKYWLPVDQYIGGIEHAIMHLMY\
FRFYHKLLRDFKLVDFDEPVKNLLCQGMVLSEAFYKIDSNSQRKWFNSSSVLIKRNTKGEIIESHTQKGEKLIYAGM\
IKMSKSKNNGIEPELIIQRYGADTIRLFIMFSAPVESDLEWKESGLKGIYRFLKKLWMLIFNYIDIKNTHKKINFDF\
LNHQQSELRYQLHKTIAKVSDDIGRRQTFNTAISEIMKLVNQLSKAPIKEEQDKSIMRESLICIIKMLYPFTPHFCF\
FVWNYFNNHSSIDNEKWPIFQKDILSKKYSTIVAQINGKKRCAIKISDSLTKEEIFLYIQNQPIIKKYLEDVDIKKI\
IYIPKKIINFVT
>YP_002468184.1
MKKNIILNLIGLRCPEPIMIIRKTIRDMKDNEKILILSDDPATKRDIPNFCYFMEHKLLKNEIKVKPYRYLLKKGL
>YP_002468032.1
MNKILVIILFSLVSLTWGTTWIAMKIATETIPPFFATGMRFLVASPILIILAYFTKTPLLFPYGQRWFQLIISIFYF\
SIPFTLMLYGGMYVSSSVASIIFSSMPVAVLTVSFLYLKKKLFLTQKIGMFISLITLLTVLLIELESQCFFQWKGIL\
ALLFALFSHAFIYSECQKKCLNVSVITFNALPSLLSGILLSTTSWFIESPHINTFSNRSILAIFYLGDFSGIFGILS\
YFYLQQKVSAFYASTVFLIFPVIAGFLENYIYKNAILLCEIWFIFPLIIGILLTLIPVNYLKKISNKITRLSCLLF
"""

seqs_prot_2 = """>YP_004590122.1
MNLEYNPKKIESFVQQYWRNNNTFSVSENSTKKKYFCVPMLPYPSGNLHMGHVRNYTISDVIARYQRMLGKNVLQPI\
GWDAFGLPAEETAIKNNISPSEWTLLNIKTMKHQLQSLGFSYDWKRELNTCNPYYYKWEQWFFIQLYKKKLVYKKNT\
LVNWCPHDKTVLANEQVHHGKCWRCNTNVILKKIPQWFMKITNYAEELLQDLKKLKNGLKKVLRMQKNWIGKSTGLL\
IKCKIYHSKKKIKIYTTKPKNIMHITFFAISMYHSIVPLLCKKNKEIHNFLNTYNNINISKMYKKNKYLGINTNHFI\
IHPVTNKKIPLWIALYIKHDYATGAIMGTPEYDKNDFHFASFHSLPINIKNLTYKNIDINEKNNSILNVHKENKDLL\
TSNDFNNTKILSNIIKILIKKKKAKKYICYKLKDWSISRQRYWGAPIPIVYTKKGKIIPVPEKNLPVILPDYSSVKN\
YTQPLQFQKSWLYTVINKENVIRETDTFDTFIESSWYYARYTCPNFTSGMIAPEKAKYWLPVDQYIGGIEHAVMHLI\
YFRFYHKLLRDFGLVSSDEPVKKLICQGMVLSDAFYTYDHNKKKIWIKPSDIKNHKKYSDKKNKNKIIYAGKIKMSK\
SKNNGIDPHIIINKYGADTLRLFLMFAAPIESALEWNDQSIIGMYRFLKKLWTFTYSLSLYLIKIPLNIHEQHDHHD\
IFIILNNTISYVTDDIKRRHSFNTAIAHIIKLFNIILKLPCQTDKNKIIIKQSLAIILKMLYPFTPHICFILWKEIY\
GKDADIDQESWPITYEIKKNIISHPFIIQINGKKKNIINISNIYSKKEKIEFALKNKKIKKYLKNKNIKNTIYIHQK\
LINFVI
>YP_004590123.1
MTELITLNLLGLRCPEPLMVLRKNIRSLKEGQIIRVLTDDFSSTRDIKIFCHFMKHILLSFSIKKIPYQYIIKIGKI\
II
>YP_004590028.1
MWGTTWIAMKIVITTIPPIFATGLRFLIIAPLSMVTAWFSDTPLLFPVGQRIMQIYISIFYFSIPFTLMLYGGRYVN\
VPTASLIFSGMPIISLLISYIIFNELINSYQFIGMSIYFISLIFFLLLQWKSSHIHQELGVLLLFISLICQSIIFIY\
FKKKFNHISVLSFNSIPSLISGILLIIFGWNIEHPVLNNFSMQSMCAICYLSVFVGFLGTISYLFLQKKIDSCYASI\
VFIIFPIVSLFLDRYMYMTKVSNFEYFFIMFLLFSVIITLFTSKKNICLIKKITKQRFNRKVA
"""

if __name__ == '__main__':
    main()



