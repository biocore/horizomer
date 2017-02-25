# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

#
# Select gene families (orthologous groups) from OrthoFinder results for
# subsequent analyses
#

import os
import sys
import click
import pandas as pd
from skbio import io


def sample_genes(ortho_groups_fp,
                 min_taxa_cutoff=10.0):
    """ Select gene families (orthologous groups) from OrthoFinder result

    Parameters
    ----------
    ortho_groups_fp : str
        orthologous groups definition file (OrthoFinder output file:
        Orthogroups.csv)
    min_taxa_cutoff : float
        minimum number (if > 1) or fraction (if <= 1) of taxa that contain
        one or more homologs of the query gene. default: 10.0

    Returns
    -------
    dict of dict of set of str
        { ogid : { taxon : set(protein(s)) }
            ogid: ID of the orthologous group (assigned by OrthoFinder)
            taxon: ID of the taxon (genome)
            protein: ID of the protein

    Raises
    ------
    ValueError
        minimum number of taxa cutoff is not a positive number

    Notes
    -----
    The purpose of this function is to pick gene families (orthologous groups,
    or OGs) identified by OrthoFinder for subsequent analyses. Specifically,
    the selected gene families will serve as input materials for phylogenetics-
    based HGT detection methods. Because phylogenetic inference requires
    adequate taxa to be statistically robust, OGs with too few members are
    discarded, and HGTs (if any) in these gene families will be considered
    unidentifiable due to inadequate taxon sampling.
    """
    ogs = pd.read_table(ortho_groups_fp, index_col=0, header=0, comment='#')
    taxa = ogs.columns
    # set minimum number of taxa cutoff
    if min_taxa_cutoff > 1:
        min_taxa_cutoff = int(min_taxa_cutoff)
    elif 0 < min_taxa_cutoff <= 1:
        min_taxa_cutoff = int(len(ogs.columns) * min_taxa_cutoff)
    else:
        raise ValueError('Error: %s is not a valid minimum taxa cutoff.'
                         % min_taxa_cutoff)
    genes = {}
    for ogid, members in ogs.iterrows():
        # drop OGs without query genome present
        if pd.notnull(members['query']):
            # drop OGs with number of taxa below cutoff
            if members.count() >= min_taxa_cutoff:
                # convert protein IDs into set
                genes[ogid] = {x: set(members[x].split(', '))
                               for x in taxa if pd.notnull(members[x])}
    return genes


def write_genes(genes,
                input_faa_dir,
                output_fa_dir,
                output_genes_fp):
    """ Write protein sequences of selected gene families to external files

    Parameters
    ----------
    genes : dict of dict of set of str
        { ogid : { taxon : set(protein(s)) }
        return value of sample_genes or filter_paralogs
    input_faa_dir : str
        directory of input protein sequences from the query genome and the
        selected taxa for comparison (FASTA format, one taxon per file)
    output_fa_dir : str
        directory to store output protein sequences (FASTA format, one gene
        family per file)
    output_genes_fp : str
        file to store a list of selected gene families and their members.
        format: gene1<tab>taxon1|protein1,taxon2|protein2,...
    """
    # generate a taxon to protein to gene family map
    #   so complicated because it is optimized for subsequent filesystem I/O
    prots = {}
    for ogid in genes:
        for taxon in genes[ogid]:
            if taxon not in prots:
                prots[taxon] = {}
            for prot in genes[ogid][taxon]:
                if prot in prots[taxon]:
                    prots[taxon][prot].add(ogid)
                else:
                    prots[taxon][prot] = set([ogid])
    # match taxa with faa filenames
    #   so complicated because OrthoFinder trims off the version number from
    #   an NCBI-style accession (e.g., GCF_012345.1 becomes GCF_012345)
    taxon2file = {}
    for fname in os.listdir(input_faa_dir):
        if fname.endswith('.faa'):
            taxon = fname.split('.')[0]
            if taxon in prots:
                taxon2file[taxon] = fname
    # read protein sequences
    seqs = {}
    for taxon in taxon2file:
        for seq in io.read(os.path.join(input_faa_dir, taxon2file[taxon]),
                           format='fasta'):
            id = seq.metadata['id']
            if id in prots[taxon]:
                seqs[id] = str(seq)
    # write protein sequences per selected gene family
    for ogid in sorted(genes):
        members = []
        with open(os.path.join(output_fa_dir, '%s.fa' % ogid), 'w') as f:
            for taxon in sorted(genes[ogid]):
                # restore taxon name from OrthoFinder-crippled form
                trutax = taxon2file[taxon][:-4]
                for prot in sorted(genes[ogid][taxon]):
                    if prot in seqs:
                        # sequence IDs are like: taxon|protein
                        f.write('>%s|%s\n%s\n' % (trutax, prot, seqs[prot]))
                        members.append('%s|%s' % (trutax, prot))
        with open(output_genes_fp, 'a') as f:
            f.write('%s\t%s\n' % (ogid, ','.join(members)))


@click.command()
@click.option('--ortho-groups-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Orthologous groups definition file (OrthoFinder output '
                   'file: Orthogroups.csv)')
@click.option('--min-taxa-cutoff', required=False, type=float, default=10.0,
              help='Minimum number (>1) or fraction (<=1) of taxa that '
                   'contain homologs of the query gene. Default: 10.0')
@click.option('--input-faa-dir', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input directory of protein sequences by taxon')
@click.option('--output-fa-dir', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Output directory of protein sequences by gene family')
@click.option('--output-genes-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output gene family list file')
def _main(ortho_groups_fp,
          min_taxa_cutoff,
          input_faa_dir,
          output_fa_dir,
          output_genes_fp):
    """ Select gene families from OrthoFinder result and write to files
    """
    genes = sample_genes(ortho_groups_fp, min_taxa_cutoff)
    if genes:
        write_genes(genes, input_faa_dir, output_fa_dir, output_genes_fp)
    sys.stdout.write('Number of gene families sampled: %s.\n' % len(genes))


if __name__ == "__main__":
    _main()
