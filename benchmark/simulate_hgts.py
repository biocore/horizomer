# ----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

#
# This script allows to simulate horizontal gene transfer by combining genes
# from raw nucleotide genomes.
#
# Input: the inputs are two genomes, donor and recipient (in GenBank format)
#        and options regarding the number of HGTs to simulate and the types of
#        HGTs allowed (orthologous replacement and novel gene acquisition
#        supported)
#
# Algorithm: all protein coding sequences (cds) for both genomes are extracted
#            using scikit-bio. For orthologous gene replacement, orthologous
#            genes are determined using OrthoFinder. A python function
#            combines genes from the two genomes (orthologous and novel).
#
# Output:    the resulting genomes in raw nucleotide and protein coding
#            sequences FASTA formats and a log file showing the nucleotide
#            positions of spiked genes
#
# Dependencies: this script requires OrthoFinder
#
# Additional information:
#   1. Emms, D. and Kelly, S. (2015). OrthoFinder: solving fundamental biases
#      in whole genome comparisons dramatically improves orthogroup inference
#      accuracy, GenomeBiology, 16:157
#

import sys
import click
from os import makedirs
from os.path import exists, basename, join, splitext, dirname, realpath
from shutil import move, rmtree
import subprocess
import time
import glob
import random

from skbio import Sequence


def extract_genbank(genbank_fp, verbose=False):
    """Extract protein coding sequences from GenBank record.

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format

    Returns
    -------
    seq: skbio.sequence.Sequence
        Sequence object
    genes: dictionary
        a dictionary of genes (CDS) and their info, with the key being the
        protein IDs and the value being a 4-element list including the
        translated sequence, the start and end positions in the genome
    """
    genes = {}
    if verbose:
        sys.stdout.write('\tParse GenBank record ...\n')
    seq = Sequence.read(genbank_fp, format='genbank')
    if verbose:
        sys.stdout.write('\t\tDone.\n')
    for feature in seq.interval_metadata._intervals:
        m = feature.metadata
        if m['type'] == 'CDS':
            protein_id = m['protein_id']
            translation = m['translation']
            strand = m['strand']
            # in scikit-bio, this number is the start location - 1
            start = feature.bounds[0][0] + 1
            end = feature.bounds[0][1]
            gene = protein_id.replace('\"', '')
            if gene not in genes:
                genes[gene] = [translation.replace(' ', '').replace('\"', ''),
                               start, end, strand]
            else:
                raise KeyError('%s already exists in dictionary' % gene)
    return seq, genes


def launch_orthofinder(proteomes_dir, output_dir, threads, verbose=False):
    """Launch OrthoFinder to report orthologous gene groups for two genomes.

    Parameters
    ----------
    proteomes_dir: string
        directory path storing FASTA coding sequences for complete genomes
    output_dir: string
        directory path storing OrthoFinder output files
    threads: integer
        number of threads to use
    verbose: boolean
        if True, run in verbose mode
    """
    if verbose:
        sys.stdout.write("\tLaunch OrthoFinder ..\n")
    orthofinder_command = ['bash', join(dirname(realpath(__file__)),
                                        'run_orthofinder.sh'),
                           '--working-dir', output_dir,
                           '--input-faa-dir', proteomes_dir,
                           '--py2-conda-env', 'wgshgt_py2',
                           '--threads', str(threads),
                           '--stdout', '/dev/null',
                           '--stderr', '/dev/null',
                           '--verbose', 'false']
    proc = subprocess.Popen(orthofinder_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print(stderr)
    if verbose:
        print(stdout)
        sys.stdout.write("\tDone\n")


def _parse_orthofinder_ids(ids_fp):
    """
    """
    ids = {}
    with open(ids_fp, 'r') as ids_f:
        for line in ids_f:
            line = line.strip().split()
            id_c = line[0].split(':')[0]
            if id_c not in ids:
                ids[id_c] = line[1]
            else:
                raise ValueError("IDs are non-unique %s" % id_c)
    return ids


def parse_orthofinder(results_dir):
    """Parse the output files of OrthoFinder for orthologous genes.

    Parameters
    ----------
    results_dir: string
        OrthoFinder results directory

    Returns
    -------
    species_ids: dictionary
        Keys are integers and values are original species FASTA filenames
    sequence_ids: dictionary
        Keys are in the form x_y (species_gene) and values are original
        accessions
    orthologous_groups: list of lists
        List of orthologous families between donor and recipient proteomes
    """
    orthologous_groups = []
    with open(
        glob.glob(
            join(
                results_dir,
                "clusters_OrthoFinder_*_id_pairs.txt"))[0], 'r') as id_pairs_f:
        # skip header lines
        for _ in range(7):
            next(id_pairs_f)
        # parse orthologous family groups
        for line in id_pairs_f:
            line = line.strip().split()
            # include only families with at least 2 orthologs
            if len(line[1:-1]) > 1:
                orthologous_groups.append(line[1:-1])
    sequence_ids = _parse_orthofinder_ids(
        join(results_dir, "SequenceIDs.txt"))
    species_ids = _parse_orthofinder_ids(
        join(results_dir, "SpeciesIDs.txt"))

    return species_ids, sequence_ids, orthologous_groups


def simulate_orthologous_rep(genes_donor,
                             seq_donor,
                             genes_recip,
                             seq_recip,
                             sequence_ids,
                             orthologous_groups,
                             orthologous_rep_prob,
                             percentage_hgts,
                             log_f):
    """Simulate orthologous replacement HGT.

    Parameters
    ----------
    genes_donor: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_donor: skbio.sequence.Sequence
        Sequence object for donor genome
    genes_recip: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_recip: skbio.sequence.Sequence
        Sequence object for recipient genome
    sequence_ids: dictionary
        Keys are in the form x_y (species_gene) and values are original
        accessions
    orthologous_groups: list of lists
        List of orthologous families between donor and recipient proteomes
    orthologous_rep_prob: float
        Probably HGT will be orthologous replacement
    percentage_hgts: float
        Percent of HGTs to simulate
    log_f: file descriptor
        Log file descriptor

    Returns
    -------
    seq_recip: skbio.sequence.Sequence
        recipient genome sequence with HGTs

    Notes
    -----
    Using list of orthologous genes between donor and recipient genomes,
    randomly choose genes to exchange from donor to recipient and output
    results to FASTA protein and nucleotide files.

    Algorithm:
        1. Choose randomly N orthogroups to be used for simulating HGTs
        2. For each orthogroup, choose randomly a donor and recipient gene
        3. Replace recipient gene with donor and output to FASTA protein
           and nucleotide files.
    """
    # number of HGTs to simulate
    num_hgts = int(percentage_hgts*orthologous_rep_prob*len(genes_recip))
    if num_hgts < 1:
        num_hgts = 1
    num_orthogroups = len(orthologous_groups)
    idx = random.sample(range(num_orthogroups), num_hgts)
    log_f.write("#type\tdonor\tstart\tend\trecipient\tnew label "
                "recipient\tstart\tend\tstrand\n")
    seq_recip_seq = str(seq_recip)
    for i in idx:
        orthogroup = orthologous_groups[i]
        substitute_genes = ['*', '*']
        # Randomly select two orthologous genes from the same family
        # representing the donor and recipient genomes. Each orthogroup is
        # guranteed to have at least two genes, one from donor (prefixed with
        # '0') and second from recipient (prefixed with '1'). The following
        # while loop will continue until an index for two genes prefixed with
        # '0' and '1' is selected. At each iteration the chances the while
        # loop must continue reduce exponentially since the index is chosen
        # randomly from the same set of options.
        while '*' in substitute_genes:
            idx2 = random.randrange(0, len(orthogroup))
            gene = orthogroup[idx2]
            if (gene.startswith('0') and substitute_genes[0] == '*'):
                substitute_genes[0] = sequence_ids[gene]
            elif (gene.startswith('1') and substitute_genes[1] == '*'):
                substitute_genes[1] = sequence_ids[gene]
        # match donor and recipient gene labels to results output by
        # OrthoFinder (in sequence_ids)
        gene_donor_label = None
        gene_recip_label = None
        if substitute_genes[0] in genes_donor:
            gene_donor_label = substitute_genes[0]
            gene_recip_label = substitute_genes[1]
        elif substitute_genes[1] in genes_donor:
            gene_donor_label = substitute_genes[1]
            gene_recip_label = substitute_genes[0]
        else:
            raise ValueError("Gene %s and %s are not in donor genome" % (
                substitute_genes[0], substitute_genes[1]))
        # rename recipient orthologous gene to donor's
        hgt_gene = "%s_hgt_o" % gene_donor_label
        genes_recip[hgt_gene] = genes_recip.pop(gene_recip_label)
        # replace recipient gene (translated sequence) with donor's
        genes_recip[hgt_gene][0] = genes_donor[gene_donor_label][0]
        # update end position of HGT gene (as it can be shorter/longer than
        # the recipient gene replaced), multiply length of substituted gene
        # by 3 to translate from codon to nucleotide length
        genes_recip[hgt_gene][2] =\
            genes_recip[hgt_gene][1] + len(genes_recip[hgt_gene][0])*3
        # replace recipient gene (nucleotide format) with donor's
        start_pos_recip, end_pos_recip, strand_recip =\
            genes_recip[hgt_gene][1:]
        start_pos_donor, end_pos_donor, strand_donor =\
            genes_donor[gene_donor_label][1:]
        seq_recip_seq = (str(seq_recip_seq[:start_pos_recip]) +
                         str(seq_donor[start_pos_donor:end_pos_donor]) +
                         str(seq_recip_seq[end_pos_recip:]))
        if strand_recip != strand_donor:
            genes_recip[hgt_gene][3] = genes_donor[gene_donor_label][3]
        # write HGTs to log file
        log_f.write("o\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            gene_donor_label,
            start_pos_donor,
            end_pos_donor,
            gene_recip_label,
            hgt_gene,
            start_pos_recip,
            end_pos_recip,
            strand_donor))
    seq_recip = Sequence(seq_recip_seq, metadata=seq_recip.metadata)
    return seq_recip


def simulate_novel_acq(genes_donor,
                       seq_donor,
                       genes_recip,
                       seq_recip,
                       orthologous_rep_prob,
                       percentage_hgts,
                       log_f):
    """Simulate novel gene acquisition HGT.

    Parameters
    ----------
    genes_donor: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_donor: skbio.sequence.Sequence
        Sequence object for donor genome
    genes_recip: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_recip: skbio.sequence.Sequence
        Sequence object for recipient genome
    orthologous_rep_prob: float
        Probably HGT will be orthologous replacement
    percentage_hgts: float
        Percent of HGTs to simulate
    log_f: file descriptor
        Log file descriptor

    Returns
    -------
    seq_recip: skbio.sequence.Sequence
        recipient genome sequence with HGTs

    Notes
    -----
    Algorithm:
        1. choose random location in recipient genome where to insert a gene
           (chosen from list of donor genes)
        2. use gene (recipient) positioning array to locate an open region
           (that doesn't include an existing gene) near the random location
           to insert the new gene (we want to avoid gene overlap so that
           compositional methods can clearly pick out individual coding
           genes)
        3. insert new gene, record existance in gene positioning array
    """
    # compute number of HGTs to simulate (novel acquisition)
    num_hgts = int(percentage_hgts*(1-orthologous_rep_prob)*len(genes_recip))
    num_hgts = max(1, num_hgts)
    # add start and end positions of recipient genome to allow for HGTs
    # simulated before the first and after the last existing gene
    gene_positions = [(0, 0), (len(seq_recip), len(seq_recip))]
    # create recipient genome gene positioning array
    for seq, start, stop, strand in genes_recip.values():
        gene_positions.append((start, stop))
    # sort array for gene positions in ascending order
    gene_positions_s = sorted(gene_positions)
    # select a random list of positions where to insert the new gene
    idx = random.sample(range(len(gene_positions_s)-1), num_hgts)
    gene_donor_labels = random.sample(list(genes_donor), num_hgts)
    log_f.write("#type\tdonor\tstart\tend\trecipient\tstart\t"
                "end\tstrand\n")
    seq_recip_seq = str(seq_recip)
    # begin simulation
    for x in range(num_hgts):
        # select donor gene (for HGT)
        gene_donor_label = gene_donor_labels[x]
        idx_recip = gene_positions_s[idx[x]][1] + 1
        # beginning from valid position for inserting new gene, check whether
        # the length can fit without overlapping with existing gene, otherwise
        # search for next valid position
        for y in range(idx[x], len(gene_positions_s)-1):
            if idx_recip + len(genes_donor[gene_donor_label][0])*3 <\
                    gene_positions_s[y+1][0]:
                # codon = sequence of 3 nucleotides (hence *3)
                idx_end = idx_recip + len(genes_donor[gene_donor_label][0])*3
                # insert gene (protein)
                hgt_gene = "%s_hgt_n" % gene_donor_label
                genes_recip[hgt_gene] =\
                    [genes_donor[gene_donor_label][0], idx_recip, idx_end,
                     genes_donor[gene_donor_label][3]]
                # insert gene (nucleotide)
                seq_recip_seq = (str(seq_recip_seq[:idx_recip]) +
                                 str(seq_donor[genes_donor[gene_donor_label]
                                               [1]:
                                               genes_donor[gene_donor_label]
                                               [2]]) +
                                 str(seq_recip_seq[idx_recip:]))
                # write HGTs to log file
                log_f.write(
                    "n\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                    (gene_donor_label,
                     genes_donor[gene_donor_label][1],
                     genes_donor[gene_donor_label][2],
                     hgt_gene,
                     idx_recip,
                     idx_end,
                     genes_donor[gene_donor_label][3]))
                break
            # try next open region
            idx_recip = gene_positions_s[y+1][1] + 1
    seq_recip = Sequence(seq_recip_seq, metadata=seq_recip.metadata)
    return seq_recip


def write_results(genes_donor,
                  donor_genome_fp,
                  genes_recip,
                  recip_genome_fp,
                  seq_donor,
                  seq_recip,
                  output_dir):
    """Write donor and recipient genomes to GenBank and FASTA files.

    Parameters
    ----------
    genes_donor: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_donor: skbio.sequence.Sequence
        Sequence object for donor genome
    genes_recip: dictionary
        A dictionary of genes, key are protein IDs values 5-element lists
    seq_recip: skbio.sequence.Sequence
        Sequence object for recipient genome

    Returns
    -------
    donor_genome_nucl_fp: string
        filepath to donor nucleotide sequence
    donor_genome_aa_fp: string
        filepath to donor protein sequences
    donor_genome_gb_fp: string
        filepath to donor genome in GenBank format
    recip_genome_nucl_fp: string
        filepath to recipient nucleotide sequence (simulated)
    recip_genome_aa_fp: string
        filepath to recipient protein sequence
    recip_genome_gb_fp: string
        filepath to donor genome in GenBank format
    """
    # output dir for simulated results
    simulated_dir = join(output_dir, "simulated")
    if not exists(simulated_dir):
        makedirs(simulated_dir)
    donor_genome_aa_fp = join(
        simulated_dir, "%s.faa" % basename(splitext(donor_genome_fp)[0]))
    recip_genome_aa_fp = join(
        simulated_dir, "%s.faa" % basename(splitext(recip_genome_fp)[0]))
    with open(donor_genome_aa_fp, 'w') as donor_genome_aa_f:
        for gene in genes_donor:
            donor_genome_aa_f.write(
                ">%s\n%s\n" % (gene, genes_donor[gene][0]))
    with open(recip_genome_aa_fp, 'w') as recip_genome_aa_f:
        for gene in genes_recip:
            recip_genome_aa_f.write(
                ">%s\n%s\n" % (gene, genes_recip[gene][0]))
    donor_genome_nucl_fp = join(
        simulated_dir, "%s.fna" % basename(splitext(donor_genome_fp)[0]))
    recip_genome_nucl_fp = join(
        simulated_dir, "%s.fna" % basename(splitext(recip_genome_fp)[0]))
    seq_donor.write(donor_genome_nucl_fp, format='fasta')
    seq_recip.write(recip_genome_nucl_fp, format='fasta')
    donor_genome_gb_fp = join(
        simulated_dir, "%s.gb" % basename(splitext(donor_genome_fp)[0]))
    recip_genome_gb_fp = join(
        simulated_dir, "%s.gb" % basename(splitext(recip_genome_fp)[0]))
    seq_donor.write(donor_genome_gb_fp, format='genbank')
    if 'LOCUS' in seq_recip.metadata:
        seq_recip.metadata['LOCUS']['size'] = len(str(seq_recip))
    seq_recip.interval_metadata._intervals = []
    for (gene, l) in sorted(genes_recip.items(), key=lambda x: x[1][1]):
        location = str(l[1]) + '..' + str(l[2])
        if l[3] == '-':
            location = 'complement(' + location + ')'
        # in scikit-bio, bounds[0][0] is the start location - 1
        bounds = [(l[1] - 1, l[2])]
        feature = {'type': 'gene', 'locus_tag': gene, '__location': location}
        seq_recip.interval_metadata.add(bounds, metadata=feature)
        feature = {'type': 'CDS', 'locus_tag': gene, '__location': location,
                   'protein_id': gene, 'translation': l[0]}
        seq_recip.interval_metadata.add(bounds, metadata=feature)
    seq_recip.write(recip_genome_gb_fp, format='genbank')
    return (donor_genome_nucl_fp, donor_genome_aa_fp, donor_genome_gb_fp,
            recip_genome_nucl_fp, recip_genome_aa_fp, recip_genome_gb_fp)


def simulate_hgts(seq_donor,
                  genes_donor,
                  seq_recip,
                  genes_recip,
                  donor_genome_fp,
                  recip_genome_fp,
                  output_dir,
                  percentage_hgts,
                  orthologous_rep_prob,
                  log_f,
                  threads=1,
                  verbose=False):
    """Simulate orthologous replacement and novel gene acquisition HGTs.

    Parameters
    ----------
    log_f: file descriptor
        Log file descriptor
    threads: integer
        number of threads to use
    """
    # output dir for OrthoFinder results
    proteomes_dir = join(output_dir, "proteomes")
    if not exists(proteomes_dir):
        makedirs(proteomes_dir)

    # write donor and recipient genes to file
    if verbose:
        sys.stdout.write("Write donor and recipient genes to file.\n")
    genes_donor_fp = join(
        proteomes_dir, "%s_donor.faa" % basename(donor_genome_fp))
    genes_recip_fp = join(
        proteomes_dir, "%s_recip.faa" % basename(recip_genome_fp))
    with open(genes_donor_fp, 'w') as genes_donor_f:
        for gene in genes_donor:
            genes_donor_f.write(">%s\n%s\n" % (gene, genes_donor[gene][0]))
    with open(genes_recip_fp, 'w') as genes_recip_f:
        for gene in genes_recip:
            genes_recip_f.write(">%s\n%s\n" % (gene, genes_recip[gene][0]))
    if verbose:
        sys.stdout.write("\tDone.\n")

    # simulate orthologous replacement
    if orthologous_rep_prob > 0.0:
        if verbose:
            sys.stdout.write("\tSimulate orthologous replacement HGTs ...\n")
        launch_orthofinder(proteomes_dir, output_dir, threads, verbose=True)
        date = time.strftime("%c").split()
        day = date[2].zfill(2)
        results_dir = join(output_dir, "orthofinder", "WorkingDirectory")
        species_ids, sequence_ids, orthologous_groups =\
            parse_orthofinder(results_dir)
        # no orthologs found, exit orthologous replacement simulation
        if orthologous_groups:
            seq_recip = simulate_orthologous_rep(genes_donor,
                                                 seq_donor,
                                                 genes_recip,
                                                 seq_recip,
                                                 sequence_ids,
                                                 orthologous_groups,
                                                 orthologous_rep_prob,
                                                 percentage_hgts,
                                                 log_f)
        else:
            print("WARNING: No orthologous genes found between donor and "
                  "recipient genome, continuing to novel gene acquisition "
                  "HGT simulation (if option selected).")
        if verbose:
            sys.stdout.write(" done.\n")
    # simulate novel gene acquisition
    if float(1-orthologous_rep_prob) > 0.0:
        if verbose:
            sys.stdout.write("\tSimulate novel gene acquisition HGTs ...")
        seq_recip = simulate_novel_acq(genes_donor,
                                       seq_donor,
                                       genes_recip,
                                       seq_recip,
                                       orthologous_rep_prob,
                                       percentage_hgts,
                                       log_f)
        if verbose:
            sys.stdout.write(" done.\n")
    return write_results(
        genes_donor=genes_donor,
        donor_genome_fp=donor_genome_fp,
        genes_recip=genes_recip,
        recip_genome_fp=recip_genome_fp,
        seq_donor=seq_donor,
        seq_recip=seq_recip,
        output_dir=output_dir)


def simulate_genbank(donor_genbank_fp,
                     recipient_genbank_fp,
                     output_dir,
                     percentage_hgts,
                     orthologous_rep_prob,
                     log_f,
                     threads,
                     verbose=False):
    """ Simulate HGTs using genuine genomes (GenBank format).

    Parameters
    ----------
    donor_genbank_fp: string
        file path to genome (donor of HGTs) in GenBank format
    recipient_genbank_fp: string
        file path to genome (recipient of HGTs) in GenBank format
    output_dir: string
        output directory path
    percentage_hgts: float
        percentage of HGT genes to simulate (of total genes in recipient
        genome)
    orthologous_rep_prob: float
        rate of orthologous replacement HGTs
    log_f: file descriptor
        Log file descriptor
    threads: integer
        number of threads to use
    """
    if verbose:
        sys.stdout.write("Parsing donor GenBank record ...\n")
    seq_donor, genes_donor = extract_genbank(donor_genbank_fp, verbose)
    if verbose:
        sys.stdout.write("\tDone.\n")
    if verbose:
        sys.stdout.write("Parsing recipient GenBank record ...\n")
    seq_recip, genes_recip = extract_genbank(recipient_genbank_fp, verbose)

    if verbose:
        sys.stdout.write("\tDone.\n")

    return simulate_hgts(
        threads=threads,
        verbose=verbose,
        seq_donor=seq_donor,
        genes_donor=genes_donor,
        seq_recip=seq_recip,
        genes_recip=genes_recip,
        donor_genome_fp=donor_genbank_fp,
        recip_genome_fp=recipient_genbank_fp,
        output_dir=output_dir,
        percentage_hgts=percentage_hgts,
        orthologous_rep_prob=orthologous_rep_prob,
        log_f=log_f)


@click.command()
@click.option('--donor-genbank-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to genome (donor of HGTs) in GenBank format')
@click.option('--recipient-genbank-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to genome (recipient of HGTs) in GenBank '
                   'format')
@click.option('--output-dir', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory path')
@click.option('--percentage-hgts', required=False, type=float, default=0.05,
              show_default=True, help='Percentage of HGT genes to simulate '
                                      '(of total genes in recipient genome)')
@click.option('--orthologous-rep-prob', required=False, type=float,
              default=0.5, show_default=True,
              help='Probability of orthologous replacement HGT (the remainder'
                   ' will be HGTs in the form of novel gene acquisition)')
@click.option('--threads', required=False, type=int, default=1,
              show_default=True, help='Number of threads to use')
@click.option('--verbose', required=False, type=bool, default=False,
              show_default=True, help='Run program in verbose mode')
def _main(donor_genbank_fp,
          recipient_genbank_fp,
          output_dir,
          percentage_hgts,
          orthologous_rep_prob,
          threads,
          verbose):
    """ Simulate HGTs by combining genes from two genomes.
    """
    if verbose:
        sys.stdout.write("Begin simulation.\n")
    if not exists(output_dir):
        makedirs(output_dir)
    log_fp = join(output_dir, "log.txt")
    with open(log_fp, 'w') as log_f:
        if verbose:
            sys.stdout.write("Simulate HGTs in genomes.\n")
        simulate_genbank(
            donor_genbank_fp=donor_genbank_fp,
            recipient_genbank_fp=recipient_genbank_fp,
            output_dir=output_dir,
            percentage_hgts=percentage_hgts,
            orthologous_rep_prob=orthologous_rep_prob,
            log_f=log_f,
            threads=threads,
            verbose=verbose)

if __name__ == "__main__":
    _main()
