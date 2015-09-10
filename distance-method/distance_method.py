#!/usr/bin/python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

""" Implement Distance Method for HGT detection based on algorithm described
    in:
        Wei. X et al., "A Distance-Based Method for Detecting HGT in Whole
        Genomes", International Symposium on Bioinformatics Research and
        Applications (ISBRA), 2008, pages 26-37

    The workflow follows the algorithm:
    1. For each gene in target genome,
        i.    BLAST sequence against all other genes in the reference genomes;
        ii.   Go to step 3 if gene has more than threshold number of homologs
              (min-num-homologs), otherwise go to next gene in target genome;
        iii.  Compute multiple sequence alignment on homolog genes using
              CLUSTAL;
        iv.   Compute pairwise distance matrix using PHYLIP's protdist
              function and Z-score normalize the set of pairwise distances
              for each gene family and species;
        v.    Add distance matrix for all pairwise distances into global
              distance matrix storing results for all genes

    2. Cluster gene families by species,
        vi.   Compute all species sets (sets of genes whose orthologs are
              detectable in exactly the same subset of the considered
              species);
        vii.  Cluster genes to each core species set cluster using the Hamming
              distance clustering algorithm;
        viii. Run outlier detection algorithm on each cluster (paragraph 2
              of section 'Detecting Outlier Genes' in original paper)
"""

import sys
import click
import numpy
import operator
import threading
import subprocess
import traceback
import shlex
from os import listdir
from os.path import join, splitext, basename
from itertools import imap


from skbio.parse.sequences import parse_fasta


class Command(object):
    """
    Enables to run subprocess commands in a different thread
    with TIMEOUT option.

    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    https://gist.github.com/kirpit/1306188
    """
    command = None
    process = None
    status = None
    output, error = '', ''

    def __init__(self, command):
        if isinstance(command, basestring):
            command = shlex.split(command)
        self.command = command

    def run(self, timeout=None, **kwargs):
        """ Run a command then return: (status, output, error). """
        def target(**kwargs):
            try:
                self.process = subprocess.Popen(self.command, **kwargs)
                self.output, self.error = self.process.communicate()
                self.status = self.process.returncode
            except:
                self.error = traceback.format_exc()
                self.status = -1
        # default stdout and stderr
        if 'stdout' not in kwargs:
            kwargs['stdout'] = subprocess.PIPE
        if 'stderr' not in kwargs:
            kwargs['stderr'] = subprocess.PIPE
        # thread
        thread = threading.Thread(target=target, kwargs=kwargs)
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.status, self.output, self.error


def hamming(str1, str2):
    "Compute the Hamming distance between two strings"
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(imap(ne, str1, str2))


def preprocess_data(working_dir,
                    target_proteomes_dir,
                    verbose=False):
    """ Map each gene to sudo name (ex. 1_1 for species 1, gene 1) for easier
    output comparison and the 10 character name limitation in PHYLIP output.
    This format is limited up to 9999 species and 99999 genes per species.

    Parameters
    ----------
    working_dir:  string
      path to working directory
    target_proteomes_dir: string
      path to directory holding proteomes for all target organisms
    verbose: boolean, optional
      output details about the running processes of this function

    Return
    ------
    gene_map: dictionary
      "two-way" dictionary storing gene names as keys and their pseudo
      names as values, and vica versa
    ref_db: dictionary
      dictionary storing FASTA label as key and sequence as value for the
      reference databases
    species: integer
      the number of species in the reference databases
    """
    gene_map = {}
    ref_db = {}
    if verbose:
        sys.stdout.write("Target organism\tNumber of genes\n")
    # each file contains genes for species
    species = 0
    for _file in listdir(target_proteomes_dir):
        if verbose:
            sys.stdout.write("%s. %s\t" % (species+1, basename(_file)))
        with open(join(target_proteomes_dir, _file), 'rb') as readfile:
            gene = 0
            for label, seq in parse_fasta(readfile):
                label = label.split()[0]
                ref_db[label] = seq
                sudo_label = "%s_%s" % (species, gene)
                gene_map[label] = sudo_label
                gene_map[sudo_label] = label
                gene += 1
        if verbose:
            sys.stdout.write("%s\n" % gene)
        species += 1
    return gene_map, ref_db, species


def launch_blast(query_proteome_fp,
                 ref_fp,
                 working_dir,
                 e_value=10e-20,
                 threads=1,
                 debug=False):
    """ Launch BLASTp given a query proteome and a reference database of
    proteomes.

    Parameters
    ----------
    query_proteome_fp: string
      filepath to query proteome
    ref_fp: string
      filepath to reference proteome
    working_dir: string
      working directory path
    e_value: float, optional
      the cutoff E-value for BLASTP results
    threads: integer
      number of threads to use for running BLASTP
    debug: boolean
      if True, run function in debug mode

    Return
    ------
    out_file_fp: string
      filepath to tabular alignment file output by
      BLASTP
    """
    db_file_fp = join(working_dir, "%s" % basename(ref_fp))
    # build blast database
    makeblastdb_command = ["makeblastdb",
                           "-in", ref_fp,
                           "-out", db_file_fp,
                           "-dbtype", "prot"]
    proc = subprocess.Popen(makeblastdb_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if (stderr and debug):
        print "[DEBUG] %s\n" % stderr

    # launch blast
    out_file_fp = join(
        working_dir, "%s.blast" % basename(splitext(query_proteome_fp)[0]))
    blastp_command = ["blastp",
                      "-db", db_file_fp,
                      "-query", query_proteome_fp,
                      "-evalue", str(e_value),
                      "-num_threads", str(threads),
                      "-outfmt", "6 std qcovs",
                      "-task", "blastp",
                      "-out", out_file_fp]
    proc = subprocess.Popen(blastp_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if (stderr and debug):
        print "[DEBUG] %s\n" % stderr

    return out_file_fp


def parse_blast(alignments_fp,
                hits,
                gene_map,
                debug=False):
    """ Parse BLASTp alignment file into a dictionary where the keys are the
    queries and the values are all the reference sequences to which the query
    mapped with E-value cutoff score.

    Parameters
    ----------
    alignments_fp: string
      filepath to tabular alignment file output by BLASTP
    hits: dictionary
      dictionary storing query (gene) names as keys and the best aligning
      reference sequences as values (one alignment per reference sequence)
    gene_map: dictionary
      "two-way" dictionary storing gene names as keys and their pseudo
      names as values, and vica versa
    debug: boolean
      if True, run function in debug mode
    """
    # read blastp results
    with open(alignments_fp, 'U') as alignments_f:
        for line in alignments_f:
            if debug:
                sys.stdout.write("[DEBUG] %s" % line)
            line = line.split()
            query = line[0]
            ref = line[1]
            if query not in hits:
                hits[query] = [ref]
            else:
                # check that the query mapped to a different species
                # since we only want the best homolog per species
                current_species = gene_map[ref].split('_')[0]
                add_alignment = True
                for gene in hits[query]:
                    species = gene_map[gene].split('_')[0]
                    if species == current_species:
                        add_alignment = False
                        break
                if add_alignment:
                    hits[query].append(ref)


def launch_msa(fasta_in_fp,
               clustal_command_fp,
               gene_map,
               ref_db,
               hits,
               query,
               timeout):
    """ Create multiple sequence alignments for all gene othologs
    using Clustalw.

    Parameters
    ----------
    fasta_in_fp: string
      filepath to FASTA file of protein sequences to use as input to
      Clustalw
    clustal_command_fp: string
      filepath to Clustalw command (interactive)
    gene_map: dictionary
      "two-way" dictionary storing gene names as keys and their pseudo
      names as values, and vica versa
    ref_db: dictionary
      dictionary storing FASTA label as key and sequence as value for the
      reference databases
    hits: dictionary
      dictionary storing query (gene) names as keys and the best aligning
      reference sequences as values (one alignment per reference sequence)
    query: string
      query gene name
    timeout: integer
      number of seconds to allow Clustalw to run before terminating the
      process
    """
    with open(fasta_in_fp, 'w') as in_f:
        sorted(hits[query])
        for ref in hits[query]:
            in_f.write(">%s\n%s\n" % (gene_map[ref], ref_db[ref]))

    with open(clustal_command_fp, 'U') as clustal_command_f:
        clustalw_command = Command("clustalw")
        status, output, error = clustalw_command.run(
            timeout=timeout,
            stdin=clustal_command_f,
            close_fds=True)
        if status < 0:
            sys.stdout.write(
                "status: %s\noutput: %s\terror: %s\t" % (
                    status, output, error))


def compute_distances(phylip_command_fp,
                      warnings=False):
    """ Compute distances between each pair of sequences in MSA (multiple
        sequence alignment) using PHYLIP's protdist function.

        Parameters
        ----------
        phylip_command_fp: string
          filepath to the PHYLIP command (interactive)
        warnings: boolean, optional
          print warnings output by PHYLIP
    """
    with open(phylip_command_fp, 'U') as phylip_command_f:
        proc = subprocess.Popen("protdist",
                                stdin=phylip_command_f,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr and warnings:
            print stderr


def normalize_distances(phylip_fp,
                        full_distance_matrix,
                        num_species,
                        full_distance_matrix_offset,
                        species_set_dict,
                        gene_bitvector_map,
                        debug=False):
    """ This function parses the output file of PHYLIP's protdist function
    containing the distance matrix, Z-score normalizes the set of pairwise
    distances between the gene in a species and all other species and stores
    the results in a separate array.

    Each normalized distance matrix is then sorted by species name and added
    to the complete array storing distance matrices for all genes. In
    addition, a list of missing species (species which did not include a
    certain gene) is also maintained and used for setting nan's in array cells
    which represent those species.

    Below is an example of a parsed distance matrix
    for 3 genes and 3 species:
        0         1         2        (genes)
    0_0 nan       nan       nan
    0_1 0.53099   0.878855  0.83673
    0_2 0.642856  1.083039  1.083039
    1_0 0.300297  0.300297  0.702003
    1_1 nan       nan       nan
    1_2 0.399722  0.379156  0.356543
    2_0 0.53099   0.53099   0.83673
    2_1 0.399722  0.399722  0.356543
    2_2 nan       nan       nan

    (species pairs)

    Example of Z-score normalized distance matrix from
    above:
        0            1           2          (genes)
    0_0 nan          nan         nan
    0_1 -1.40548346  0.83861735  0.56686611
    0_2 -1.41421356  0.70710678  0.70710678
    1_0 -0.70710678 -0.70710678  1.41421356
    1_1 nan          nan         nan
    1_2 1.20493966   0.03869341 -1.24363308
    2_0 -0.70710678 -0.70710678  1.41421356
    2_1 0.70710678   0.70710678 -1.41421356
    2_2 nan          nan         nan

    (species pairs)

    Parameters
    ----------
    phylip_fp: string
        filepath to distance matrix output by PHYLIP's protdist function
    full_distance_matrix: dictionary
        complete distance matrix for pairwise alignments between all species
        for every gene
    num_species: integer
        number of species in the reference database
    full_distance_matrix_offset: integer
        the index offset for elements in full_distance_matrix where to write
        the next array
    species_set_dict: dictionary
        dictionary containing the binary indicator vectors as keys and the
        number of genes with identical species set represented by the binary
        vectors as values
    gene_bitvector_map: list
        list containing the binary indicator vector for each query gene
    debug: boolean
        if True, run function in debug mode
    """
    # assume a pairwise alignment exists for all species
    missing_species = [str(x) for x in range(0, num_species)]
    # scan through file and remove species that exist
    # from missing_species list
    with open(phylip_fp, 'U') as phylip_f:
        for line in phylip_f:
            alignment_dist = line.strip().split()
            if len(alignment_dist) == 1:
                continue
            if not line.startswith(' '):
                species = line.split()[0].split('_')[0]
                missing_species.remove(species)

    # scan through file again, collecting alignment
    # distances
    orig_order_labels = []
    p = numpy.empty(shape=(num_species, num_species))
    p.fill(numpy.nan)
    idx_p = 0
    ind_a = 0
    with open(phylip_fp, 'U') as phylip_f:
        alignment_list = []
        # skip first line containing number of lines in
        # the file
        next(phylip_f)
        for line in phylip_f:
            if debug:
                sys.stdout.write("[DEBUG] %s" % line)
            alignment_dist = line.strip().split()
            if line.startswith(' '):
                alignment_list.extend(alignment_dist)
            else:
                # new species alignment pairs
                if alignment_list:
                    for i in range(0, len(missing_species)):
                        alignment_list.append(None)
                    a = numpy.asarray(alignment_list[1:], dtype=float)
                    a[ind_a] = numpy.nan
                    ind_a += 1
                    mean = numpy.nanmean(a)
                    stdev = numpy.nanstd(a)
                    # modify array with normalized values
                    for x in numpy.nditer(a, op_flags=['readwrite']):
                        x[...] = (x-mean)/stdev
                    p[idx_p] = a
                    idx_p += 1
                    orig_order_labels.append(alignment_list[0])
                alignment_list = alignment_dist

    # add distance on final line
    for i in range(0, len(missing_species)):
        alignment_list.append(None)
    a = numpy.asarray(alignment_list[1:], dtype=float)
    a[ind_a] = numpy.nan
    mean = numpy.nanmean(a)
    stdev = numpy.nanstd(a)
    for x in numpy.nditer(a, op_flags=['readwrite']):
        x[...] = (x-mean)/stdev
    p[idx_p] = a
    orig_order_labels.append(alignment_list[0])

    # add the missing species names to the labels array
    bitvector_gene = 'I' * num_species
    for species in missing_species:
        orig_order_labels.append("%s_X" % species)
        # indicate missing gene for current species
        l = list(bitvector_gene)
        l[int(species)] = 'O'
        bitvector_gene = ''.join(l)

    # update species set counts
    if bitvector_gene not in species_set_dict:
        species_set_dict[bitvector_gene] = 1
    else:
        species_set_dict[bitvector_gene] += 1

    gene_bitvector_map[full_distance_matrix_offset] = bitvector_gene

    # sort the distance matrix based on species names (S1, S2, S3 ..)
    # in order to be consistent across all gene families
    sorted_order_labels = sorted(orig_order_labels)
    map_orig_sorted = {}
    for idx_orig, label in enumerate(orig_order_labels):
        idx_sorted = sorted_order_labels.index(label)
        map_orig_sorted[idx_orig] = idx_sorted

    # re-order rows and columns by ordered species name (0,1,2 ..)
    p2 = numpy.zeros(shape=(num_species, num_species))
    for idx_a, arr in enumerate(p):
        t = numpy.zeros(shape=num_species)
        for idx_b, el in enumerate(arr):
            t[map_orig_sorted[idx_b]] = el
        p2[map_orig_sorted[idx_a]] = t
    del p

    # add normalized distance matrix for current gene
    # to full distance matrix
    full_distance_matrix[full_distance_matrix_offset] = p2


def cluster_distances(species_set_dict,
                      species_set_size,
                      hamming_distance):
    """ Cluster gene families by species with detectable orthologs in exactly
    the same subset of the considered species.

    Ex. Assume we have 4 genes and 5 species with the following distance
    matrix:

        0             1             2             3
    0_0 nan           nan           nan           nan
    0_1 -1.59564844   -1.388031632  -0.9634704748 -1.342272936
    0_2 -0.4259542606 nan           0.7035215923  1.223837777
    0_3 -1.55041393   -1.51499567   -0.9634704748 -1.330178178
    0_4 -0.3659762821 0.8346037464  0.6565725705  1.15274682
    ..
    ..

    There are two binary indicator vectors to represent the species present in
    the four genes: IIIII (gene 0, 2 and 3), II0II (gene 1). If the core set
    threshold was 3, then there would be 1 core species set represented by
    IIIII.

    Parameters
    ----------
    species_set_dict: dictionary
      dictionary containing the binary indicator vectors as
      keys and the number of genes with identical species
      set represented by the binary vectors as values
    species_set_size: integer
      threshold number of genes in a species set to
      allow it to form a core cluster
    hamming_distance: integer
      maximum number of mismatches between two binary
      indicator vectors (ex. IIII and I0II) for the
      genes in a candidate vector to be merged into the
      core cluster

    Returns
    -------
    gene_clusters_dict: dictionary
      dictionary containing core species sets as keys
      and all belonging species sets as values
      (determined by the Hamming distance clustering
      algorithm)
    """
    sorted_species_set = sorted(species_set_dict.items(),
                                key=operator.itemgetter(1), reverse=True)

    # determine core clusters (initial species sets
    # with more than species_set_size genes)
    gene_clusters_dict = {}
    # if the largest species set contains less than
    # threshold (species_set_size) elements, set the
    # only core cluster to the largest species set
    if sorted_species_set[0][1] < species_set_size:
        gene_clusters_dict[sorted_species_set[0][0]] = []
    else:
        for bitvector in sorted_species_set:
            if bitvector[1] >= species_set_size:
                gene_clusters_dict[bitvector[0]] = []

    # assign species sets with fewer than species_set_size
    # species to core clusters if the Hamming distance
    # between the two bitvectors is less than
    # hamming_distance
    species_set_assigned = []
    for cluster_core in gene_clusters_dict:
        for bitvector in sorted_species_set:
            if (bitvector[0] not in species_set_assigned and
                    hamming(cluster_core, bitvector[0]) <= hamming_distance):
                gene_clusters_dict[cluster_core].append(bitvector[0])
                species_set_assigned.append(bitvector[0])

    # assign the remaining species sets to the
    # cluster with the closest core Hamming distance
    for bitvector in sorted_species_set:
        if bitvector[0] not in species_set_assigned:
            min_hamming_cluster = ""
            min_hamming_distance = sys.maxint
            # find cluster core with smallest Hamming distance
            # to species set
            for cluster_core in gene_clusters_dict:
                dist = hamming(cluster_core, bitvector[0])
                if dist < min_hamming_distance:
                    min_hamming_distance = dist
                    min_hamming_cluster = cluster_core
            gene_clusters_dict[min_hamming_cluster].append(bitvector[0])

    return gene_clusters_dict


def detect_outlier_genes(species_set,
                         gene_bitvector_map,
                         full_distance_matrix,
                         stdev_offset,
                         outlier_hgt,
                         num_species,
                         total_genes,
                         debug=False):
    """ Detect outlier genes algorithm described in section "Detecting
    `Outlier' Genes" of the Wei. X et al. paper. The full distance matrix is
    represented in the format:

    full_distance_matrix[#genes][#species][#species] =
    [[[0_0, 0_1, 0_2, .., 0_n]
      [1_0, 1_1, 1_2, .., 1_n]
      ..
      [n_0, n_1, n_2, .., n_n]]

     [[0_0, 0_1, 0_2, .., 0_n]
      [1_0, 1_1, 1_2, .., 1_n]
      ..
      [n_0, n_1, n_2, .., n_n]]

      ..
     [[0_0, 0_1, 0_2, .., 0_n]
      [1_0, 1_1, 1_2, .., 1_n]
      ..
      [n_0, n_1, n_2, .., n_n]]]

    The mean and standard deviation are computed for each species pair
    including all genes.

    Parameters
    ----------
    species_set: list
        list of bitvectors representing species clusters to use in detecting
        outlier genes
    gene_bitvector_map: list
        list containing the binary indicator vector for each query gene
    full_distance_matrix: dictionary
        complete distance matrix for pairwise alignments between all species
        for every gene
    stdev_offset: integer
        the number of standard deviations a gene's normalized distance is from
        the mean to identify it as an outlier for a species pair
    outlier_hgt: float
        the fraction (value between (0,1]) of normalized pairwise distances
        over all species-pair vectors belonging to the same gene that are
        z-score standard deviations from the mean
    num_species: integer
        number of species in the reference database
    total_genes: integer
        total number of genes in the query genome with at least
        min_num_homologs (determined by BLAST search)
    debug: boolean
        if True, run function in debug mode

    Returns
    -------
    outlier_genes: set
        set of atypical genes
    """
    numpy.around(full_distance_matrix, decimals=5, out=full_distance_matrix)
    outlier_flag_matrix = numpy.zeros(
        shape=(total_genes, num_species, num_species), dtype=bool)
    distance_vector = numpy.zeros(total_genes)
    if debug:
        sys.stdout.write("[DEBUG] species_species\t")
        for k in range(total_genes):
            sys.stdout.write("gene # %s".ljust(12) % k)
        sys.stdout.write("[low_bound, up_bound]\n")
    for i in range(num_species):
        for j in range(num_species):
            if i != j:
                for k in range(total_genes):
                    distance_vector[k] = full_distance_matrix[k][i][j]
                mean = numpy.nanmean(distance_vector)
                stdev = numpy.nanstd(distance_vector)
                low_bound = round(mean - stdev_offset*stdev, 5)
                up_bound = round(mean + stdev_offset*stdev, 5)
                if debug:
                    sys.stdout.write("[DEBUG] %s_%s\t".ljust(20) % (i, j))
                for k, distance in enumerate(distance_vector):
                    spaces = "".ljust(2)
                    if distance < 0:
                        spaces = "".ljust(1)
                    if (distance != numpy.nan and
                       ((distance < low_bound) or (distance > up_bound))):
                        outlier_flag_matrix[k][i][j] = 1
                        if debug:
                            sys.stdout.write(
                                "%s\033[92m%s\033[0m" % (spaces, distance))
                    elif debug:
                        sys.stdout.write("%s%s" % (spaces, distance))
                if debug:
                    sys.stdout.write("\t[%s, %s]\n" % (low_bound, up_bound))

    # traverse outlier_matrix by gene and count the number of
    # outlier distances by species
    outlier_count_matrix = numpy.zeros(
        shape=(total_genes, num_species), dtype=int)
    for i in range(total_genes):
        for j in range(num_species):
            for k in range(num_species):
                if outlier_flag_matrix[i][j][k]:
                    outlier_count_matrix[i][k] += 1

    # if number of outlier distances exceeds threshold, label
    # gene as outlier
    outlier_genes = set()
    for i in range(total_genes):
        for j in range(num_species):
            if outlier_count_matrix[i][j] > num_species*outlier_hgt:
                outlier_genes.add(i)

    return outlier_genes


def output_full_matrix(matrix, num_species):
    """ Output distance matrix to stdout
    """
    for i in range(num_species):
        for j in range(num_species):
            # for gene number
            for k in range(len(matrix)):
                sys.stdout.write("%s\t" % matrix[k][i][j])
            sys.stdout.write("\n")


@click.command()
@click.argument('query-proteome-fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('target-proteomes-dir', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('working-dir', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.option('--min-num-homologs', type=int, required=False, default=3,
              show_default=True, help="The mininum number of homologs "
                                      "(determined by BLAST search) for each "
                                      "gene to test")
@click.option('--e-value', type=float, required=False, default=10e-20,
              show_default=True, help="The E-value cutoff to identify "
                                      "orthologous genes using BLASTP")
@click.option('--threads', type=int, required=False, default=1,
              show_default=True, help="Number of threads to use")
@click.option('--stdev-offset', type=float, required=False, default=2.326,
              show_default=True, help="The number of standard deviations a "
                                      "gene's normalized distance is from "
                                      "the mean to identify it as an outlier "
                                      "for a species pair")
@click.option('--outlier-hgt', type=float, default=0.5, show_default=True,
              required=False, help="The fraction (value between (0,1]) of "
                                   "normalized pairwise distances over all "
                                   "species-pair vectors belonging to the "
                                   "same gene that are z-score standard "
                                   "deviations from the mean")
@click.option('--species-set-size', type=int, required=False, default=30,
              show_default=True, help="Threshold number of genes to consider "
                                      "a species set large (a species set is "
                                      "a set of genes whose orthologs are "
                                      "detectable in exactly the same subset "
                                      "of the considered species)")
@click.option('--hamming-distance', type=int, required=False, default=2,
              show_default=True, help="Distance between two binary vectors "
                                      "indicating the species in which the "
                                      "corresponding ortholog gene appears")
@click.option('--verbose', type=bool, required=False, default=False,
              show_default=True, help="Run in verbose mode")
@click.option('--debug', type=bool, required=False, default=False,
              show_default=True, help="Run in debug mode")
@click.option('--warnings', type=bool, required=False, default=False,
              show_default=True, help="Print program warnings")
@click.option('--timeout', type=int, required=False, default=120,
              show_default=True, help="Number of seconds to allow Clustalw "
                                      "to run per call")
def distance_method_main(query_proteome_fp,
                         target_proteomes_dir,
                         working_dir,
                         min_num_homologs,
                         e_value,
                         threads,
                         stdev_offset,
                         outlier_hgt,
                         species_set_size,
                         hamming_distance,
                         verbose,
                         debug,
                         warnings,
                         timeout):
    """ Main function to run the Distance-Method HGT detection algorithm

    Parameters
    ----------
    query_proteome_fp: string
        file path to query proteome (assumed HGT recipient)
    target_proteomes_dir: string
        dir path to target proteomes (assumed HGT donors)
    working_dir: string
        dir path to working directory
    min_num_homologs: integer
        mininum number of homologs for each gene to test
    e_value: float
        E-value cutoff to identify orthologous genes using BLASTP
    threads: integer
        number of threads to use
    stdev_offset: float
        number of standard deviations a gene's normalized distance is from the
        mean to identify it as an outlier for a species pair
    outlier_hgt: float
        the fraction (value between (0,1]) of normalized pairwise distances
        over all species-pair vectors belonging to the same gene that are
        Z-score standard deviations from the mean
    species_set_size: integer
        threshold number of genes to consider a species set large
    hamming_distance: integer
        distance between two binary vectors indicating the species in which
        the corresponding ortholog gene appears
    verbose: boolean
        run in verbose mode
    debug: boolean
        run in debug mode
    warnings: boolean
        output warnings
    timeout: integer
        number of seconds to allow Clustalw to run per call
    """
    if verbose:
        sys.stdout.write(
            "Begin whole-genome HGT detection using the Distance method.\n\n")
        sys.stdout.write("Query genome: %s\n" % query_proteome_fp)

    gene_map, ref_db, num_species = preprocess_data(
        working_dir=working_dir,
        target_proteomes_dir=target_proteomes_dir,
        verbose=verbose)

    if debug:
        sys.stdout.write("\n[DEBUG] gene map:\n")
        for gene in gene_map:
            sys.stdout.write("[DEBUG] %s: %s\n" % (gene, gene_map[gene]))

    if verbose:
        sys.stdout.write("\nRunning BLASTp ..\n")
    hits = {}

    for _file in listdir(target_proteomes_dir):
        # launch BLASTp
        alignments_fp = launch_blast(
            query_proteome_fp=query_proteome_fp,
            ref_fp=join(target_proteomes_dir, _file),
            working_dir=working_dir,
            e_value=e_value,
            threads=threads,
            debug=debug)

        # generate a dictionary of orthologous genes
        parse_blast(alignments_fp=alignments_fp,
                    hits=hits,
                    gene_map=gene_map,
                    debug=debug)

    # keep only genes with >= min_num_homologs
    hits_min_num_homologs = {}
    max_homologs = 0
    for query in hits:
        len_hits = len(hits[query])
        if query in hits[query]:
            len_hits -= 1
        if len_hits >= min_num_homologs:
            if query in hits_min_num_homologs:
                raise ValueError("Duplicate gene names found: %s" % query)
            hits_min_num_homologs[query] = hits[query]
            if len_hits > max_homologs:
                max_homologs = len_hits
    hits.clear()

    if verbose:
        sys.stdout.write(
            "Total number of orthologous gene families with at "
            "least %s genes: %s\n" % (
                min_num_homologs, len(hits_min_num_homologs)))
    if debug:
        sys.stdout.write("[DEBUG] Blast matches:\n")
        for query in hits_min_num_homologs:
            sys.stdout.write(
                "[DEBUG] %s: %s\n" % (query, hits_min_num_homologs[query]))
    # generate command for CLUSTALW
    phy_msa_fp = join(working_dir, "msa.phy")
    dnd_msa_fp = join(working_dir, "msa.dnd")
    phylip_fp = join(working_dir, "msa.dis")
    # create fasta file for each gene family and run CLUSTALW
    fasta_in_fp = join(working_dir, "input.faa")
    clustal_command_fp = join(working_dir, "clustal_command.txt")
    with open(clustal_command_fp, 'w') as clustal_command_f:
        clustal_command_f.write(
            '1\n%s\n2\n9\n1\n4\n\n1\n%s\n%s\nX\n\nX\n' % (
                fasta_in_fp, phy_msa_fp, dnd_msa_fp))
    phylip_command_fp = join(working_dir, "phylip_command.txt")
    with open(phylip_command_fp, 'w') as phylip_command_f:
        phylip_command_f.write('%s\nF\n%s\nR\nY\n' % (phy_msa_fp, phylip_fp))

    total_genes = len(hits_min_num_homologs)
    if verbose:
        sys.stdout.write("\nRunning CLUSTALW and PROTDIST ..\n")
    if max_homologs > num_species:
        raise ValueError(
            "max_homologs > num_species: %s > %s " % (
                max_homologs, num_species))
    # distance matrix containing distances between all ortholog genes
    full_distance_matrix = numpy.zeros(
        shape=(total_genes, num_species, num_species), dtype=float)
    # dictionary to store all subsets of orthologs (keys) and
    # their number of occurrences (values) (maximum occurrences
    # is equal to the number of genes)
    species_set_dict = {}
    gene_bitvector_map = {}
    gene_id = {}
    for i, query in enumerate(hits_min_num_homologs):
        if verbose:
            print "Computing MSA and distances for gene %s .. (%s/%s)" % (
                query, i+1, total_genes)
        gene_id[i] = query
        # generate a multiple sequence alignment
        # for each orthologous gene family
        launch_msa(fasta_in_fp=fasta_in_fp,
                   clustal_command_fp=clustal_command_fp,
                   ref_db=ref_db,
                   gene_map=gene_map,
                   hits=hits_min_num_homologs,
                   query=query,
                   timeout=timeout)

        # compute distances between each pair of sequences in MSA
        compute_distances(phylip_command_fp=phylip_command_fp,
                          warnings=warnings)

        # Z-score normalize distance matrix and add results
        # to full distance matrix (for all genes)
        normalize_distances(phylip_fp=phylip_fp,
                            full_distance_matrix=full_distance_matrix,
                            num_species=num_species,
                            full_distance_matrix_offset=i,
                            species_set_dict=species_set_dict,
                            gene_bitvector_map=gene_bitvector_map,
                            debug=debug)

    # output_full_matrix(full_distance_matrix, num_species)

    # cluster gene families by species
    gene_clusters_dict = cluster_distances(
        species_set_dict=species_set_dict,
        species_set_size=species_set_size,
        hamming_distance=hamming_distance)

    # detect outlier genes per core cluster of genes
    sys.stdout.write("\nCandidate HGT genes: \n")
    for core_cluster in gene_clusters_dict:
        outlier_genes = detect_outlier_genes(
            species_set=gene_clusters_dict[core_cluster],
            gene_bitvector_map=gene_bitvector_map,
            full_distance_matrix=full_distance_matrix,
            stdev_offset=stdev_offset,
            outlier_hgt=outlier_hgt,
            num_species=num_species,
            total_genes=total_genes,
            debug=debug)

        for gene in outlier_genes:
            sys.stdout.write("%s\n" % gene_id[gene])

    # output_full_matrix(outlier_genes, num_species)


if __name__ == "__main__":
    distance_method_main()
