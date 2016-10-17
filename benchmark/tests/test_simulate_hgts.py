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
from os import makedirs
from os.path import join, dirname, abspath
import time
import copy
from operator import itemgetter

from skbio import Sequence
import skbio.io

from benchmark.simulate_hgts import (extract_genbank,
                                     launch_orthofinder,
                                     parse_orthofinder,
                                     simulate_orthologous_rep,
                                     simulate_novel_acq,
                                     write_results,
                                     simulate_genbank)


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

        self.genes_donor = {'D_1': ['MEKEYSPKKIENYVQEFWKK', 30, 90, '+'],
                            'D_2': ['MKKNIILNLIGLRCPEPIMI', 140, 200, '+'],
                            'D_3': ['MNKILVIILFSLVSLTWGTTWIAMKIA',
                                    220, 301, '-']}
        self.seq_donor = Sequence(
            "TTACAGAATTTGTAGATCCATGTATTTTTGATGGAGAAGGAGTACAGCCCCAA"
            "GAAGATCGAGAACTACGTGCAGGAGTTCTGGAAGAAGTTTATATATCTTTTTT"
            "TATTAATAATTAATATAAAAAAATTTTCTATATTATGAAGAAGAACATCATCC"
            "TGAACCTGATCGGCCTGAGGTGCCCCGAGCCCATCATGATCCGTAGTCGTGAT"
            "GCTGATGTATGAACAAGATCCTGGTGATCATCCTGTTCAGCCTGGTGAGCCTG"
            "ACCTGGGGCACCACCTGGATCGCCATGAAGATCGCCTGTATATACCAGCAATT"
            "TCCCCAATTTTGCTTTTAAAATTTGAAATTGATTTTTTTATTTTAGAAAACGT"
            "TGGTTTTTGACCAGTAATATATTTTATTGA")
        self.seq_donor.metadata['LOCUS'] = {'locus_name': 'donor',
                                            'size': len(str(self.seq_donor)),
                                            'unit': 'bp',
                                            'shape': 'circular',
                                            'division': 'CON',
                                            'mol_type': 'DNA',
                                            'date': '01-JAN-1900'}
        self.genes_recip = {'R_1': ['MNLEYNPKKIESFVQQYWRN', 45, 104, '+'],
                            'R_2': ['MTELITLNLLGLRCPEPLMV', 120, 179, '+'],
                            'R_3': ['MWGTTWIAMKIVITTIPPIFATGLRFL',
                                    240, 320, '+']}
        self.seq_recip = Sequence(
            "AATTAGCCTATTAAATTTATTAATTTTTATATTGCTTAATATATAATGAATTT"
            "GGAATATAATCCAAAAAAAATTGAATCTTTTGTTCAACAATATTGGAGAAATT"
            "TAATATTAATTTGAATGACTGAATTGATTACTTTGAATTTGTTGGGTTTGAGA"
            "TGTCCAGAACCATTGATGGTTTTTTTTTACTAATTTTTATTATTAACAAAAAA"
            "TTTAAAGATATTTTTAAAAAATTCAATAATGTGGGGTACTACTTGGATTGCTA"
            "TGAAAATTGTTATTACTACTATTCCACCAATTTTTGCTACTGGTTTGAGATTT"
            "TTGATTTTAGACTGGTATTTCAAGAACGATTACTTTTAAACTGGCGTTTAAAT"
            "ATCAACATCTCCCAGCTATCCTACACAAAAA")
        self.seq_recip.metadata['LOCUS'] = {'locus_name': 'recipient',
                                            'size': len(str(self.seq_recip)),
                                            'unit': 'bp',
                                            'shape': 'circular',
                                            'division': 'CON',
                                            'mol_type': 'DNA',
                                            'date': '01-JAN-1900'}

    def tearDown(self):
        rmtree(self.working_dir)

    def test_extract_genbank(self):
        """Test parsing sequence and gene information from a GenBank record.
        """
        seq, genes = extract_genbank(
            genbank_fp=join(self.root, "genbank_sample_record.gbk"),
            verbose=True)
        genes_exp = {'AAA98665.1': ['SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKR'
                                    'AVVSSASEAAEVLLRVDNIIRARPRTANRQHM', 0,
                                    206, '+'],
                     'AAA98667.1': ['MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYT'
                                    'TYQSFNLPQFVPINRHPALIDYIEELILDVLSKLT'
                                    'HVYRFSICIINKKNDLCIEKYVLDFSELQHVDKDD'
                                    'QIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTIT'
                                    'FEAVINAIELELGHKLDRNRRVDSLEEKAEIERDS'
                                    'NWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPL'
                                    'IIHQFSEKLISGDDKILNGVYSQYEEGESIFGSLF',
                                    3299, 4037, '-'],
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
                                    686, 3158, '+']}
        self.assertEqual(str(seq), sample_seq)
        self.assertDictEqual(genes, genes_exp)

    def test_launch_orthofinder(self):
        """Test running OrthoFinder.
        """
        launch_orthofinder(self.proteomes_dir, 1, verbose=True)
        date = time.strftime("%c").split()
        day = date[2]
        if int(day) < 10:
            day = "0%s" % day
        results_dir = join(
            self.proteomes_dir, "Results_%s%s" % (date[1], day))
        orthogroups_exp = [['YP_002468181.1', 'YP_004590122.1'],
                           ['YP_004590123.1', 'YP_002468184.1'],
                           ['YP_002468032.1', 'YP_004590028.1']]
        orthogroups_act = []
        with open(join(results_dir, "OrthologousGroups.txt"), 'r') as o:
            orthogroups_act = [line.split()[1:] for line in o]
        orthogroups_act_sorted = [sorted(group) for group in orthogroups_act]
        orthogroups_exp_sorted = [sorted(group) for group in orthogroups_exp]
        orthogroups_act_sorted_2 = sorted(orthogroups_act_sorted)
        orthogroups_exp_sorted_2 = sorted(orthogroups_exp_sorted)
        self.assertListEqual(orthogroups_act_sorted_2,
                             orthogroups_exp_sorted_2)

    def test_parse_orthofinder(self):
        """Test parsing OrthoFinder results.
        """
        launch_orthofinder(self.proteomes_dir, 1)
        date = time.strftime("%c").split()
        day = date[2]
        if int(day) < 10:
            day = "0%s" % day
        results_dir = join(
            self.proteomes_dir, "Results_%s%s" % (date[1], day),
            "WorkingDirectory")
        species_ids, sequence_ids, orthogroups_act =\
            parse_orthofinder(results_dir)
        species_ids_exp = {'1': 'seqs_prot_2.fasta', '0': 'seqs_prot_1.fasta'}
        seq_ids_exp = {'1_2': 'YP_004590028.1', '1_1': 'YP_004590123.1',
                       '1_0': 'YP_004590122.1', '0_2': 'YP_002468032.1',
                       '0_0': 'YP_002468181.1', '0_1': 'YP_002468184.1'}
        orthogroups_exp = [['0_0', '1_0'], ['0_1', '1_1'], ['0_2', '1_2']]
        self.assertDictEqual(species_ids, species_ids_exp)
        self.assertDictEqual(sequence_ids, seq_ids_exp)
        orthogroups_act_sorted = [sorted(group) for group in orthogroups_act]
        orthogroups_exp_sorted = [sorted(group) for group in orthogroups_exp]
        orthogroups_act_sorted_2 = sorted(orthogroups_act_sorted)
        orthogroups_exp_sorted_2 = sorted(orthogroups_exp_sorted)
        self.assertListEqual(orthogroups_act_sorted_2,
                             orthogroups_exp_sorted_2)

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
        with open(log_fp, 'r') as log_f:
            hgts_sim = [line.strip().split()
                        for line in log_f
                        if not line.startswith('#')]
        # verify HGTs
        for hgt in hgts_sim:
            recip_label_replaced = hgt[4]
            self.assertTrue(recip_label_replaced not in self.genes_recip)
            hgt_g = hgt[5]
            self.assertTrue(hgt_g in self.genes_recip)
            donor_label = hgt[1]
            # protein sequences are equal for donor and recipient gene
            self.assertEqual(genes_donor_orig[donor_label][0],
                             self.genes_recip[hgt_g][0])
            recip_start = hgt[6]
            recip_end = hgt[7]
            # length of gene (nucleotide format) is consistent with the
            # gene start and end positions on the recipient genome
            self.assertEqual(len(self.genes_recip[hgt_g][0])*3,
                             int(recip_end)-int(recip_start))
            # genome subsequence representing HGT in recipient genome is
            # equal to the nucleotide format of the donor gene
            self.assertEqual(
                str(seq_recip[
                    self.genes_recip[hgt_g][1]:self.genes_recip[hgt_g][2]]),
                translated_nucl[donor_label])

    def test_simulate_novel_acq(self):
        """Test simulating novel gene acquisition HGTs.
        """
        # genes_recip and seq_recip are modified within function
        # simulate_orthologous_rep, store their original copies for
        # downstream verification
        genes_donor_orig = copy.deepcopy(self.genes_donor)
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
        self.assertEqual(len(self.genes_donor), len(genes_donor_orig))
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
        with open(log_fp, 'r') as log_f:
            hgts_sim = [line.strip().split()
                        for line in log_f
                        if not line.startswith('#')]
        gene_positions =\
            [(self.genes_recip[g][1], self.genes_recip[g][2], True)
             if 'hgt' in g
             else (self.genes_recip[g][1], self.genes_recip[g][2], False)
             for g in self.genes_recip]
        # verify HGTs
        for hgt in hgts_sim:
            donor_label = hgt[1]
            hgt_l = hgt[4]
            recip_start = int(hgt[5])
            recip_end = int(hgt[6])
            self.assertTrue(hgt_l in self.genes_recip)
            self.assertEqual(recip_start, self.genes_recip[hgt_l][1])
            self.assertEqual(recip_end, self.genes_recip[hgt_l][2])
            self.assertEqual(hgt[7], self.genes_recip[hgt_l][3])
            # length of gene (nucleotide format) is consistent with the
            # gene start and end positions on the recipient genome
            self.assertEqual(len(self.genes_recip[hgt_l][0])*3,
                             int(recip_end)-int(recip_start))
            # genome subsequence representing HGT in recipient genome is
            # equal to the nucleotide format of the donor gene
            self.assertEqual(
                str(seq_recip[
                    self.genes_recip[hgt_l][1]:self.genes_recip[hgt_l][2]]),
                translated_nucl[donor_label])
        gene_positions_s = sorted(gene_positions, key=itemgetter(0))
        # verify none of the original genes overlap with HGTs
        for x in range(0, len(gene_positions_s)):
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
        dnr_g_nucl_fp, dnr_g_aa_fp, dnr_g_gb_fp, \
        rcp_g_nucl_fp, rcp_g_aa_fp, rcp_g_gb_fp =\
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
        locus = {'unit': 'bp', 'shape': 'circular', 'division': 'CON',
                 'mol_type': 'DNA', 'date': '01-JAN-1900'}
        donor_gb = Sequence.read(dnr_g_gb_fp, format='genbank')
        locus['locus_name'] = 'donor'
        locus['size'] = len(str(self.seq_donor))
        # test for correctness of donor GenBank file
        self.assertEqual(str(donor_gb), str(self.seq_donor))
        self.assertDictEqual(donor_gb.metadata['LOCUS'], locus)
        recip_gb = Sequence.read(rcp_g_gb_fp, format='genbank')
        locus['locus_name'] = 'recipient'
        locus['size'] = len(str(self.seq_recip))
        # test for correctness of recipient GenBank file
        self.assertEqual(str(recip_gb), str(self.seq_recip))
        self.assertDictEqual(recip_gb.metadata['LOCUS'], locus)
        donor_aa_dict = {}
        for seq in skbio.io.read(dnr_g_aa_fp, format='fasta'):
            donor_aa_dict[seq.metadata['id']] = str(seq)
        # test for correctness of donor protein coding sequences
        self.assertTrue(len(donor_aa_dict), len(self.genes_donor))
        for gene in donor_aa_dict:
            self.assertTrue(gene in self.genes_donor)
            self.assertEqual(donor_aa_dict[gene], self.genes_donor[gene][0])
        recip_aa_dict = {}
        for seq in skbio.io.read(rcp_g_aa_fp, format='fasta'):
            recip_aa_dict[seq.metadata['id']] = str(seq)
        # test for correctness of recipient protein coding sequences
        self.assertTrue(len(recip_aa_dict), len(self.genes_recip))
        for gene in recip_aa_dict:
            self.assertTrue(gene in self.genes_recip)
            self.assertEqual(recip_aa_dict[gene], self.genes_recip[gene][0])

    def load_seqs(self, file_fp):
        """Load FASTA file into dictionary
        """
        gene_dict = {}
        for seq in skbio.io.read(file_fp, format='fasta'):
            seq_id = seq.metadata['id']
            if seq_id not in gene_dict:
                gene_dict[seq_id] = seq
            else:
                raise ValueError("Duplicate gene %s" % seq_id)
        return gene_dict

    def test_simulate_genbank(self):
        """Test simulating HGTs using input donor and recipient GenBank files.
        """
        donor_genbank_fp = join(self.root, "genbank",
                                "GCF_000010365.1_ASM1036v1_genomic.gbff")
        recip_genbank_fp = join(self.root, "genbank",
                                "GCF_000441575.1_ASM44157v1_genomic.gbff")
        output_dir = self.simulated_dir
        percentage_hgts = 0.05
        orthologous_rep_prob = 0.5
        log_fp = join(self.working_dir, "log.txt")
        threads = 1
        with open(log_fp, 'w') as log_f:
            dnr_nucl_fp, dnr_aa_fp, dnr_gb_fp, \
            rcp_nucl_fp, rcp_aa_fp, rcp_gb_fp =\
                simulate_genbank(donor_genbank_fp,
                                 recip_genbank_fp,
                                 output_dir,
                                 percentage_hgts,
                                 orthologous_rep_prob,
                                 log_f,
                                 threads,
                                 verbose=True)
        # load all simulated HGT information from log file
        hgts_sim_ortho = {}
        hgts_sim_replc = {}
        with open(log_fp, 'r') as log_f:
            for line in log_f:
                if not line.startswith('#'):
                    line = line.strip().split()
                    hgt_type = line[0]
                    donor = line[1]
                    if hgt_type == 'o':
                        hgts_sim = hgts_sim_ortho
                    elif hgt_type == 'n':
                        hgts_sim = hgts_sim_replc
                    else:
                        raise ValueError(
                            "HGT type %s not supported" % hgt_type)
                    if donor not in hgts_sim:
                            hgts_sim[donor] = line[2:]
                    else:
                        raise ValueError(
                            "Duplicate gene donor %s" % donor)
        # load sequences from donor proteome
        dnr_prot_genes = self.load_seqs(dnr_aa_fp)
        # load sequences from simulated recipient proteome
        rcp_prot_hgts = self.load_seqs(rcp_aa_fp)
        # check HGTs in recipient proteome match sequences from the donor,
        # and if HGT was by orthologous replacement, the replaced gene does
        # not exist in recipient proteome
        num_simulated_hgts = 0
        for seq_id_recip in rcp_prot_hgts:
            if "hgt" in seq_id_recip:
                num_simulated_hgts += 1
                seq_id_donor = seq_id_recip.split('_hgt_')[0]
                self.assertEqual(
                    str(rcp_prot_hgts[seq_id_recip]),
                    str(dnr_prot_genes[seq_id_donor]))
                if "hgt_o" in seq_id_recip:
                    replaced_gene = hgts_sim_ortho[seq_id_donor][2]
                    self.assertTrue(replaced_gene not in rcp_prot_hgts)


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
