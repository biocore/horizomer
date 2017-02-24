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

import sys
import click
import pandas as pd


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

    Return
    ------
    dict of dict of set of str
        { ogid : { taxon : set(protein(s)) }
            ogid: ID of the orthologous group (assigned by OrthoFinder)
            taxon: ID of the taxon (genome)
            protein: ID of the protein

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


@click.command()
@click.option('--ortho-groups-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Orthologous groups definition file (OrthoFinder output '
                   'file: Orthogroups.csv)')
@click.option('--min-taxa-cutoff', required=False, type=float, default=10.0,
              help='Minimum number (>1) or fraction (<=1) of taxa that '
                   'contain homologs of the query gene. Default: 10.0')
def _main(ortho_groups_fp,
          min_taxa_cutoff):
    """ Select gene families from OrthoFinder result and write to files
    """
    genes = sample_genes(ortho_groups_fp, min_taxa_cutoff)
    sys.stdout.write(str(genes))


if __name__ == "__main__":
    _main()
