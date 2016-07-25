# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

#
# Parse output files of GI tools.
#

import sys
import click
from skbio import Sequence


def parse_output_gi(genbank_fp,
                    gi_fp):
    """ Identify genes contained in GIs

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format
    gi_fp: string
        file path to GI list

    Returns
    -------
    output: string
        gene names (protein_ids) separated by newline
    """
    genes = {}
    gb = Sequence.read(genbank_fp, format='genbank')
    for feature in gb.interval_metadata.features:
        if feature['type_'] == 'CDS':
            if 'protein_id' in feature:
                protein_id = feature['protein_id'].replace('\"', '')
                loc = gb.interval_metadata.features[feature]
                start = loc[0][0] + 1
                end = loc[0][1]
                if protein_id not in genes:
                    genes[protein_id] = [start, end]
    genes_in_gi = {}
    with open(gi_fp, 'r') as input_f:
        for line in input_f:
            l = line.strip().split()
            if len(l) < 2:
                continue
            (start, end) = (int(l[0]), int(l[1]))
            for (gene, l) in genes.items():
                if (l[0] >= start and l[1] <= end):
                    if gene not in genes_in_gi:
                        genes_in_gi[gene] = 1
    output = []
    for gene in sorted(genes_in_gi):
        output.append(gene)
    return '\n'.join(output)


@click.command()
@click.option('--genbank-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Genome in GenBank format')
@click.option('--gi-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory path')
def _main(genbank_fp,
          gi_fp):
    """ Identify genes contained in GIs
    """
    output = parse_output_gi(genbank_fp, gi_fp)
    sys.stdout.write(output)


if __name__ == "__main__":
    _main()
