#!/usr/bin/python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Implement Distance Method for HGT detection based on algorithm described
#    in:
#        Wei. X et al., "A Distance-Based Method for Detecting HGT in Whole
#        Genomes", International Symposium on Bioinformatics Research and
#        Applications (ISBRA), 2008, pages 26-37
#
#    The workflow follows the algorithm:
#    1. For each gene in target genome,
#        i.    BLAST sequence against all other genes in the reference
#              genomes;
#        ii.   Go to step 3 if gene has more than threshold number of homologs
#              (min-num-homologs), otherwise go to next gene in target genome;
#        iii.  Compute multiple sequence alignment on homolog genes using
#              CLUSTAL;
#        iv.   Compute pairwise distance matrix using PHYLIP's protdist
#              function and Z-score normalize the set of pairwise distances
#              for each gene family and species;
#        v.    Add distance matrix for all pairwise distances into global
#              distance matrix storing results for all genes
#
#    2. Cluster gene families by species,
#        vi.   Compute all species sets (sets of genes whose orthologs are
#              detectable in exactly the same subset of the considered
#              species);
#        vii.  Cluster genes to each core species set cluster using the
#              Hamming distance clustering algorithm;
#        viii. Run outlier detection algorithm on each cluster (paragraph 2
#              of section 'Detecting Outlier Genes' in original paper)
#
#    Requires protdist version 3.696
#

import sys
import click
import numpy
import operator
import threading
import subprocess
import traceback
import shlex
from os.path import join, basename, isdir, exists, getsize
from os import mkdir

from glob import glob

import skbio.io


class Command(object):
    """Run subprocess commands in a different thread with TIMEOUT option.

    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    https://gist.github.com/kirpit/1306188
    """
    process = None
    status = None
    output, error = '', ''

    def __init__(self, command):
        if isinstance(command, str):
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
    """Compute the Hamming distance between two strings.

    Parameters
    ----------
    str1: string
        string
    str2: string
        string
    """
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))


def preprocess_data(working_dir,
                    target_proteomes_dir,
                    extensions,
                    verbose=False):
    """ Map each gene to sudo name (ex. 1_1 for species 1, gene 1).

    Parameters
    ----------
    working_dir:  string
        path to working directory
    target_proteomes_dir: string
        path to directory holding proteomes for all target organisms
    extensions: list
        list of extensions for reference proteomes
    verbose: boolean, optional
        output details about the running processes of this function

    Returns
    -------
    gene_map: dictionary
        "two-way" dictionary storing gene names as keys and their pseudo
        names as values, and vica versa
    ref_db: dictionary
        dictionary storing FASTA label as key and sequence as value for the
        reference databases
    species: integer
        the number of species in the reference databases

    Notes
    -----
        This will facilitate easier output comparison and the 10 character
        name limitation in PHYLIP output. This format is limited up to 9999
        species and 99999 genes per species.
    """
    gene_map = {}
    ref_db = {}
    if verbose:
        sys.stdout.write("Target organism\tNumber of genes\n")
    # each file contains genes for species
    files = [f
             for ext in extensions
             for f in glob("%s/*%s" % (target_proteomes_dir, ext))]
    for species, _file in enumerate(files):
        if verbose:
            sys.stdout.write("%s. %s\t" % (
                species+1, basename(_file)))
        for gene, seq in enumerate(skbio.io.read(_file, format='fasta')):
            label = seq.metadata['id']
            ref_db[label] = seq
            sudo_label = "%s_%s" % (species, gene)
            if label in gene_map:
                raise ValueError("Duplicate sequence labels are "
                                 "not allowed: %s" % label)
            gene_map[label] = sudo_label
            gene_map[sudo_label] = label
        if verbose:
            sys.stdout.write("%s\n" % gene)
    return gene_map, ref_db, species+1


def launch_diamond(query_proteome_fp,
                   ref_fp,
                   working_dir,
                   tmp_dir,
                   e_value=10e-20,
                   threads=1,
                   debug=False):
    """ Launch DIAMOND for a query and a reference database of proteomes.

    Parameters
    ----------
    query_proteome_fp: string
      filepath to query proteome
    ref_fp: string
      filepath to reference proteome
    working_dir: string
      working directory path
    tmp_dir:
      temporary working directory for DIAMOND
    e_value: float, optional
      the cutoff E-value for BLASTP results
    threads: integer
      number of threads to use for running DIAMOND BLASTP
    debug: boolean
      if True, run function in debug mode

    Returns
    -------
    out_file_fp: string
      filepath to tabular alignment file output by DIAMOND
    """
    db_file_fp = join(working_dir, "%s" % basename(ref_fp))
    # build DIAMOND database
    makediamonddb_command = ["diamond",
                             "makedb",
                             "--in", ref_fp,
                             "-d", db_file_fp,
                             "--threads", str(threads)]
    proc = subprocess.Popen(makediamonddb_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if (stderr and debug):
        print("[DEBUG] %s\n" % stderr)

    # launch DIAMOND
    out_file_fp = join(
        working_dir, "%s.daa" % basename(query_proteome_fp))
    diamond_command = ["diamond",
                       "blastp",
                       "-t", tmp_dir,
                       "--db", "%s.dmnd" % db_file_fp,
                       "--query", query_proteome_fp,
                       "--evalue", str(e_value),
                       "--threads", str(threads),
                       "--daa", out_file_fp,
                       "--sensitive"]
    proc = subprocess.Popen(diamond_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if (stderr and debug):
        print("[DEBUG] %s\n" % stderr)

    # convert output to tab delimited file
    out_file_conv_fp = join(
        working_dir, "%s.m8" % basename(query_proteome_fp))
    diamond_convert_command = ["diamond",
                               "view",
                               "--daa", out_file_fp,
                               "-f", "tab",
                               "-o", out_file_conv_fp]
    proc = subprocess.Popen(diamond_convert_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if (stderr and debug):
        print("[DEBUG] %s\n" % stderr)

    return out_file_conv_fp


def launch_blast(query_proteome_fp,
                 ref_fp,
                 working_dir,
                 e_value=10e-20,
                 threads=1,
                 debug=False):
    """ Launch BLASTp for a query and a reference database of proteomes.

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

    Returns
    -------
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
        print("[DEBUG] %s\n" % stderr)

    # launch blast
    out_file_fp = join(
        working_dir, "%s.blast" % basename(query_proteome_fp))
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
        print("[DEBUG] %s\n" % stderr)

    return out_file_fp


def parse_blast(alignments_fp,
                hits,
                gene_map,
                debug=False):
    """ Parse BLASTp alignment file into a dictionary.

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

    Notes
    -----
        The keys are the queries and the values are all the reference
        sequences to which the query mapped with E-value cutoff score.
    """
    # read blastp results
    with open(alignments_fp, 'r') as alignments_f:
        for line in alignments_f:
            if debug:
                sys.stdout.write("[DEBUG] %s" % line)
            query, ref = line.split()[:2]
            if query not in hits:
                hits[query] = [ref]
            else:
                # check that the query mapped to a different species
                # since we only want the best homolog per species
                if gene_map[ref].split('_')[0] not in [
                        gene_map[gene].split('_')[0] for gene in hits[query]]:
                    hits[query].append(ref)


def launch_msa(fasta_in_fp,
               clustal_command_fp,
               gene_map,
               ref_db,
               hits,
               query,
               timeout):
    """ Create MSA for all gene othologs using Clustalw.

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
        for ref in hits[query]:
            in_f.write(">%s\n%s\n" % (gene_map[ref], ref_db[ref]))

    with open(clustal_command_fp, 'r') as clustal_command_f:
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
    """ Compute distances between each pair of sequences in the MSA.

    Parameters
    ----------
    phylip_command_fp: string
      filepath to the PHYLIP command (interactive)
    warnings: boolean, optional
      print warnings output by PHYLIP

    Notes
    -----
        Use PHYLIP's protdist function.
    """
    with open(phylip_command_fp, 'r') as phylip_command_f:
        proc = subprocess.Popen("protdist",
                                stdin=phylip_command_f,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr and warnings:
            print(stderr)


def normalize_distances(phylip_fp,
                        full_distance_matrix,
                        num_species,
                        full_distance_matrix_offset,
                        species_set_dict,
                        gene_bitvector_map,
                        debug=False):
    """ Parse and normalize the output file of PHYLIP's protdist function.

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

    Notes
    -----
        Parse PHYLIP's protdist output containing the distance matrix, Z-score
        normalize the set of pairwise distances between the gene in a species
        and all other species and stores the results in a separate array.

        Each normalized distance matrix is then sorted by species name and
        added to the complete array storing distance matrices for all genes.
        In addition, a list of missing species (species which did not include
        a certain gene) is also maintained and used for setting nan's in array
        cells which represent those species.

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
    """
    # assume a pairwise alignment exists for all species
    missing_species = [str(x) for x in range(0, num_species)]
    # scan through file and remove species that exist
    # from missing_species list
    if exists(phylip_fp) and getsize(phylip_fp) > 0:
        with open(phylip_fp, 'r') as phylip_f:
            next(phylip_f)
            for line in phylip_f:
                if not line.startswith(' '):
                    species = line.split()[0].split('_')[0]
                    missing_species.remove(species)
    else:
        raise ValueError('%s does not exist or is empty' % phylip_fp)

    # scan through file again, collecting alignment
    # distances
    orig_order_labels = []
    p = numpy.empty(shape=(num_species, num_species))
    p.fill(numpy.nan)
    idx = 0
    with open(phylip_fp, 'r') as phylip_f:
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
                    a[idx] = numpy.nan
                    p[idx] = (a - numpy.nanmean(a)) / numpy.nanstd(a)
                    idx += 1
                    orig_order_labels.append(alignment_list[0])
                alignment_list = alignment_dist

    # add distance on final line
    for i in range(0, len(missing_species)):
        alignment_list.append(None)
    a = numpy.asarray(alignment_list[1:], dtype=float)
    a[idx] = numpy.nan
    p[idx] = (a - numpy.nanmean(a)) / numpy.nanstd(a)
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
    x = numpy.argsort(numpy.array(orig_order_labels))

    # re-order rows and columns by ordered species name (0,1,2 ..)
    p2 = numpy.zeros(shape=(num_species, num_species))
    for idx_a, arr in enumerate(p):
        t = numpy.zeros(shape=num_species)
        for idx_b, el in enumerate(arr):
            t[x[idx_b]] = el
        p2[x[idx_a]] = t
    del p

    # add normalized distance matrix for current gene
    # to full distance matrix
    full_distance_matrix[full_distance_matrix_offset] = p2


def cluster_distances(species_set_dict,
                      species_set_size,
                      hamming_distance):
    """ Hamming distance clustering algorithm

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
    gene_clusters_list: list of tuples
        list of tuples containing core species sets and all belonging species
        sets (determined by the Hamming distance clustering algorithm)

    Notes
    -----
        Cluster gene families by species with detectable orthologs in exactly
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

        There are two binary indicator vectors to represent the species
        present in the four genes: IIIII (gene 0, 2 and 3), II0II (gene 1). If
        the core set threshold was 3, then there would be 1 core species set
        represented by IIIII.
    """
    sorted_species_set = sorted(list(species_set_dict.items()),
                                key=operator.itemgetter(1), reverse=True)
    # determine core clusters (initial species sets with more than
    # species_set_size genes)
    gene_clusters_list = []
    # if the largest species set contains less than threshold
    # (species_set_size) elements, set the only core cluster to the largest
    # species set
    if sorted_species_set[0][1] < species_set_size:
        cluster_core = (sorted_species_set[0][0], [])
        gene_clusters_list = []
    for bitvector in sorted_species_set:
        if bitvector[1] >= species_set_size:
            cluster_core = (bitvector[0], [])
            gene_clusters_list.append(cluster_core)
    # assign species sets with fewer than species_set_size species to core
    # clusters if the Hamming distance between the two bitvectors is less than
    # hamming_distance
    species_set_assigned = []
    for idx, cluster_core in enumerate(gene_clusters_list):
        for bitvector in sorted_species_set:
            bv = bitvector[0]
            if (bv not in species_set_assigned and
                    hamming(cluster_core[0], bv) <= hamming_distance):
                gene_clusters_list[idx][1].append(bv)
                species_set_assigned.append(bv)
    # assign the remaining species sets to the cluster with the closest core
    # Hamming distance
    for bitvector in sorted_species_set:
        bv = bitvector[0]
        if bv not in species_set_assigned:
            min_hamming_cluster = -1
            min_hamming_distance = sys.maxsize
            # find cluster core with smallest Hamming distance to species set
            for idx, cluster_core in enumerate(gene_clusters_list):
                dist = hamming(cluster_core[0], bv)
                if dist < min_hamming_distance:
                    min_hamming_distance = dist
                    min_hamming_cluster = idx
            if min_hamming_cluster >= 0:
                gene_clusters_list[min_hamming_cluster][1].append(bv)

    return gene_clusters_list


def detect_outlier_genes(species_set,
                         gene_bitvector_map,
                         full_distance_matrix,
                         stdev_offset,
                         outlier_hgt,
                         num_species,
                         total_genes,
                         debug=False):
    """ Detect outlier genes.

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

    Notes
    -----
        Algorithm described in section "Detecting `Outlier' Genes" of the Wei.
        X et al. paper. The full distance matrix is represented in the format:

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

    # traverse outlier_matrix by gene and count the number of outlier
    # distances by species
    outlier_count_matrix = numpy.zeros(
        shape=(total_genes, num_species), dtype=int)
    for i in range(total_genes):
        for j in range(num_species):
            for k in range(num_species):
                if outlier_flag_matrix[i][j][k]:
                    outlier_count_matrix[i][k] += 1

    # if number of outlier distances exceeds threshold, label gene as outlier
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


def distance_method(query_proteome_fp,
                    target_proteomes_dir,
                    working_dir,
                    output_hgt_fp,
                    align_software,
                    tabular_alignments_fp=None,
                    ext=['fa', 'fasta', 'faa'],
                    min_num_homologs=3,
                    e_value=10e-20,
                    threads=1,
                    stdev_offset=2.326,
                    outlier_hgt=0.5,
                    species_set_size=30,
                    hamming_distance=2,
                    verbose=False,
                    debug=False,
                    warnings=False,
                    timeout=120):
    """ Run Distance Method algorithm

    Parameters
    ----------
    query_proteome_fp: string
        filepath to query proteome
    target_proteomes_dir: string
        dirpath to target proteomes
    working_dir: string
        dirpath to working directory
    output_hgt_fp: string
        filepath to output file for storing detected HGTs
    align_software: string
        software to use for sequence alignment (BLAST or DIAMOND)
    tabular_alignments_fp: string, optional
        filepath to tabular sequence alignments
    ext: list, optional
        list of file extensions to open in the target proteomes directory
    min_num_homologs: integer, optional
        the mininum number of homologs (determined by BLAST search)
        for each gene to test
    e_value: float, optional
        the E-value cutoff to identify orthologous genes using BLASTP
    threads: integer, optional
        number of threads to use for sequence alignment
    stdev_offset: float, optional
        the number of standard deviations a gene's normalized distance
        is from the mean to identify it as an outlier for a species pair
    outlier_hgt: float, optional
        the fraction (value between (0,1]) of normalized pairwise distances
        over all species-pair vectors belonging to the same gene that are
        z-score standard deviations from the mean
    species_set_size: integer, optional
        threshold number of genes to consider a species set large (a species
        set is a set of genes whose orthologs are detectable in exactly the
        same subset of the considered species)
    hamming_distance: integer, optional
        distance between two binary vectors indicating the species in which
        the corresponding ortholog gene appears
    verbose: boolean, optional
        if True, run in verbose mode
    debug: boolean, optional
        if True, run in debug mode
    warnings: boolean, optional
        if True, output warnings
    timeout: integer, optional
        number of seconds to allow Clustalw to run per call
    """
    if verbose:
        sys.stdout.write(
            "Begin whole-genome HGT detection using the Distance method.\n\n")
        sys.stdout.write("Query genome: %s\n" % query_proteome_fp)

    extensions = set(['fa', 'fasta', 'faa'])
    extensions.update(ext)

    # create working directory if doesn't exist
    if not isdir(working_dir):
        mkdir(working_dir)

    gene_map, ref_db, num_species = preprocess_data(
        working_dir=working_dir,
        target_proteomes_dir=target_proteomes_dir,
        extensions=extensions,
        verbose=verbose)

    if debug:
        sys.stdout.write("\n[DEBUG] gene map:\n")
        for gene in gene_map:
            sys.stdout.write("[DEBUG] %s: %s\n" % (gene, gene_map[gene]))

    if verbose:
        sys.stdout.write("\nRunning BLASTp ..\n")
    hits = {}

    # tabular alignments provided
    if tabular_alignments_fp is not None:
        # generate a dictionary of orthologous genes
        parse_blast(alignments_fp=tabular_alignments_fp,
                    hits=hits,
                    gene_map=gene_map,
                    debug=debug)
    # tabular alignments to be created
    else:
        files = [f
                 for e in extensions
                 for f in glob("%s/*%s" % (target_proteomes_dir, e))]
        for _file in files:
            # launch BLASTp
            if align_software == "blast":
                alignments_fp = launch_blast(
                    query_proteome_fp=query_proteome_fp,
                    ref_fp=_file,
                    working_dir=working_dir,
                    e_value=e_value,
                    threads=threads,
                    debug=debug)
            elif align_software == "diamond":
                alignments_fp = launch_diamond(
                    query_proteome_fp=query_proteome_fp,
                    ref_fp=_file,
                    working_dir=working_dir,
                    tmp_dir=working_dir,
                    e_value=e_value,
                    threads=threads,
                    debug=debug)
            else:
                raise ValueError(
                    "Software not supported: %s" % align_software)

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
    open(phy_msa_fp, 'a').close()
    dnd_msa_fp = join(working_dir, "msa.dnd")
    open(dnd_msa_fp, 'a').close()
    phylip_fp = join(working_dir, "msa.dis")
    open(phylip_fp, 'a').close()
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
            print("Computing MSA and distances for gene %s .. (%s/%s)" % (
                query, i+1, total_genes))
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
    with open(output_hgt_fp, 'w') as output_hgt_f:
        output_hgt_f.write("\n# Candidate HGT genes: \n")
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

            if outlier_genes:
                for gene in outlier_genes:
                    output_hgt_f("%s\n" % gene_id[gene])

    # output_full_matrix(outlier_genes, num_species)


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
@click.argument('output-hgt-fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.option('--align-software', type=click.Choice(['diamond', 'blast']),
              required=False, default=['diamond'], show_default=True,
              help="Software to use for blasting sequences")
@click.option('--tabular-alignments-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help="Tabular alignments in m6 format (output from BLAST or "
                   "DIAMOND)")
@click.option('--ext', multiple=True, type=str, required=False,
              default=['fa', 'fasta', 'faa'], show_default=True,
              help="File extensions of target proteomes (multiple extensions "
                   "can be given by calling --ext ext1 --ext ext2)")
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
                         output_hgt_fp,
                         align_software,
                         tabular_alignments_fp,
                         ext,
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
    """ Run the Distance-Method HGT detection algorithm.
    """
    distance_method(query_proteome_fp=query_proteome_fp,
                    target_proteomes_dir=target_proteomes_dir,
                    working_dir=working_dir,
                    output_hgt_fp=output_hgt_fp,
                    align_software=align_software,
                    tabular_alignments_fp=tabular_alignments_fp,
                    ext=ext,
                    min_num_homologs=min_num_homologs,
                    e_value=e_value,
                    threads=threads,
                    stdev_offset=stdev_offset,
                    outlier_hgt=outlier_hgt,
                    species_set_size=species_set_size,
                    hamming_distance=hamming_distance,
                    verbose=verbose,
                    debug=debug,
                    warnings=warnings,
                    timeout=timeout)


if __name__ == "__main__":
    distance_method_main()
