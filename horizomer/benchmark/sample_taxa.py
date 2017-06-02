#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

#
# Select taxa (genomes) from BLASTP/DIAMOND hit table for subsequent analyses
#

import click
import pandas as pd


def sample_taxa(hit_table_fp,
                prot2tax_dict_fp):
    """ Select taxa (genomes) from BLASTP/DIAMOND hit table

    Parameters
    ----------
    hit_table_fp : str
        sequence similarity search hit table in standard tabular format
    prot2tax_dict_fp : str
        protein to taxon (genome) dictionary file

    Return
    ------
    set of str
        IDs of sampled taxa (genomes)

    Notes
    -----
    The basic idea of this script is to pick reference genomes which contains
    one or more protein-coding genes which shares sequence similarity with
    the query genome. These genomes, along with the query genome, will be
    passed to the downstream phylogenetics pipeline (OrthoFinder, Phylomizer
    and PhyloPhlAn) and the phylogenetics-based HGT detection tools.

    In order to convincingly describe the evolutionary history of an HGT event,
    sufficient taxa from both the donor group and the recipient group need to
    be included.

    In the following example, the phylogenetic tree of gene family A is
    displayed, by which a horizontal transfer of gene A from some member of
    genus 1 (likely species 1.3) to species 2.2 (***) of genus 2 is evidently
    suggested.

                                      /---species 1.1
                                 /---|
                            /---|     \---species 1.2
                           |
                           |            /-species 2.2 ***
                  /-----genus 1      /-|
                 |         |     /--|   \-species 1.3
                 |         |    |   |
                 |          \---|    \----species 1.4
                 |              |
        ---------|              |     /---species 1.5
                 |               \---|
                 |                    \---species 1.6
                 |
                 |          /-------------species 2.1
                 |         |
                  \-----genus 2       /---species 2.3
                           |     /---|
                            \---|     \---species 2.4
                                |
                                 \--------species 2.5

    The current status is a simple prototype. It indiscriminatingly samples all
    taxa indicated by subject proteins in the DIAMOND hit table (i.e., all
    genomes sharing more or less some homology with the query genome).

    In the future this code may require heavy and repeated improvement based on
    statistical and empirical studies of benchmarking results.
    """
    # column names in standard BLAST tabular format
    m8cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_table(hit_table_fp, index_col=None, names=m8cols, comment='#')
    # get all subject protein IDs mentioned in the hit table
    prots = set(df['sseqid'])
    # get IDs of taxa hosting subject proteins, based on a dictionary
    taxa = set()
    with open(prot2tax_dict_fp, 'r') as f:
        for line in f:
            # dictionary format: protein ID <tab> taxon IDs separated by comma
            prot, taxon_list = line.rstrip().split('\t')
            if prot in prots:
                taxa.update(taxon_list.split(','))
    return taxa


@click.command()
@click.option('--hit-table-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Sequence similarity search hit table in standard tabular '
                   'format (e.g., -outfmt 6 for BLAST or m8 for DIAMOND).')
@click.option('--prot2tax-dict-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Protein to taxon (genome) dictionary file.')
@click.option('--output-taxa-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output list of selected taxa.')
def _main(hit_table_fp,
          prot2tax_dict_fp,
          output_taxa_fp):
    """ Select taxa for subsequent analyses
    """
    taxa = sample_taxa(hit_table_fp, prot2tax_dict_fp)
    # write selected taxon IDs
    with open(output_taxa_fp, 'w') as f:
        for taxon in sorted(taxa):
            f.write('%s\n' % taxon)


if __name__ == "__main__":
    _main()
