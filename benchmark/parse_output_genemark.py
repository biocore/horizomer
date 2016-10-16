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


def parse_output_genemark(genbank_fp,
                          genemark_output_fp):
    """ Extract atypical genes identified by GeneMark

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format
    genemark_output_fp: string
        file path to GeneMark output gene list (*.lst)

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
                if protein_id not in genes:
                    strand = '-' if feature['rc_'] else '+'
                    loc = gb.interval_metadata.features[feature]
                    start = loc[0][0] + 1
                    end = loc[0][1]
                    genes[protein_id] = (start, end, strand)
    atypical_genes = []
    with open(genemark_output_fp, 'r') as input_f:
        reading = False
        for line in input_f:
            l = line.strip().split()
            if len(l) == 2 and l == ['#', 'Length']:
                reading = True;
            elif reading and len(l) == 6 and l[5] == '2':
                (start, end, strand) = (int(l[2].lstrip('<>')),
                                        int(l[3].lstrip('<>')),
                                        l[1])
                for (gene, l) in genes.items():
                    if l[0] == start and l[1] == end and l[2] == strand:
                        atypical_genes.append(gene)
    return '\n'.join(sorted(atypical_genes))


@click.command()
@click.option('--genbank-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Genome in GenBank format')
@click.option('--genemark-output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='GeneMark output file (*.lst)')
def _main(genbank_fp,
          genemark_output_fp):
    """ Extract atypical genes identified by GeneMark
    """
    output = parse_output_genemark(genbank_fp, genemark_output_fp)
    sys.stdout.write(output)


if __name__ == "__main__":
    _main()
