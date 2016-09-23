# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

#
# Parse output files of HGT tools.
#

import click
import sys


# T-REX version 3.6
# RANGER-DTL-U version 1.0
# RIATA-HGT version 3.5.6
# JANE version 4
# each tuple consists of three strings, first string is the unique string to
# identify the line with HGT information, second and third strings are the
# bounds for the actual number of HGTs
hgt_parse_strs = {
    'ranger-dtl': ('The minimum reconciliation cost is: ',
                   'Transfers: ',
                   ', Losses'),
    'trex': ('hgt : number of HGT(s) found = ',
             'hgt : number of HGT(s) found = ',
             ' '),
    'jane4': ('Host Switch: ',
              'Host Switch: ',
              ' '),
    'riata-hgt': ('There are ',
                  'There are ',
                  ' component(s)')}


def parse_hgts(input_f, method):
    """ Extract number of HGTs found.

    Parameters
    ----------
    input_f: string
        file descriptor for T-REX output results
    method: string
        HGT detection method

    Returns
    -------
    number_of_hgts: string
        number of HGTs reported by a tool, or NaN if an entry was not found
    """
    for line in input_f:
        if hgt_parse_strs[method][0] in line:
            return line.strip().split(
                hgt_parse_strs[method][1])[1].split(
                    hgt_parse_strs[method][2])[0]
    return 'NaN'


def parse_consel(input_f):
    """ Parse output of Consel version 0.20.

    Parameters
    ----------
    input_f: string
        file descriptor for Consel output results

    Returns
    -------
    pvalues: list
        list of P-values
    """
    pvalues = []
    # skip header lines
    skip_lines = 3
    for s in range(skip_lines):
        next(input_f)

    for line in input_f:
        line = line.split()
        # skip empty line at bottom of file
        if not line:
            continue
        pv_au = line[4]
        if 0 <= float(pv_au) <= 1:
            pvalues.append("%.2f" % float(pv_au))
    return pvalues


def parse_darkhorse(input_f, output_fp, low_lpi, high_lpi):
    """ Parse output of DarkHorse (smry file).

    Paramters
    ---------
    input_f: string
        file descriptor for Consel output results
    low_lpi: float
        lower LPI (lineage probability index) score bound
    high_lpi: float
        upper LPI score bound
    output_fp: str
        Filepath to output best hit genome IDs

    Returns
    -------
    hgts: string
        one putative HGT-derived gene per line
        columns: query_id, besthit_id, tax_id, species, lineage, pct_id,
        pct_coverage, norm_LPI

    Notes
    -----
    Parse output of DarkHorse to return tab-separated file of putative HGTs
    using the LPI bounds and a file with all best hit genome IDs.
    """
    best_hit_ids = set()
    hgts = []
    # skip header
    next(input_f)
    for l in input_f:
        l = line.strip('\r\n').split('\t')
        bets_hit_ids.add(l[3])
        if (l[5] > low_lpi) and
                (l[5] < high_lpi):
            hgt = '\t'.join(l[0], l[3], l[12], l[13], l[14], l[6], l[9], l[4])
            hgts.append(hgt)
    if output_fp:
        with open(output_fp, 'w') as output_f:
            output_f.write('\n'.join(best_hit_ids)) 
    return '\n'.join(hgts)


def parse_hgtector(input_f):
    """ Parse output of HGTector version 0.2.1.

    Parameters
    ----------
    input_f: string
        HGTector working directory path

    Returns
    -------
    output: string
        one putative HGT-derived gene per line
        columns: query_id, donor_taxid, donor_species, donor_lineage, pct_id,
        pct_coverage
    """
    hgts = []
    for line in input_f:
        l = line.strip('\r\n').split('\t')
        if (len(l) == 15) and (l[7] == '1'):
            hgt = '\t'.join((l[0], l[12], l[13], l[14], l[10], l[11]))
            hgts.append(hgt)
    return '\n'.join(hgts)


def parse_output(hgt_results_fp,
                 method,
                 low_lpi,
                 high_lpi,
                 output_fp=None):
    """Call parse_hgts() based on HGT detection method used.

    Parameters
    ----------
    hgt_results_fp: str
        filepath to detected HGTs
    method: str
        tool used to detect HGTs
    output_fp: str
        output file storing best hit IDs (DarkHorse)
    low_lpi: float
        lower bound LPI score (DarkHorse Lineage Probability Index)
    high_lpi: float
        upper bound LPI score (DarkHorse Lineage Probability Index)

    Returns
    -------
    output: string
        number of HGTs detected
    """
    with open(hgt_results_fp, 'r') as input_f:
        if (method == 'ranger-dtl' or
                method == 'trex' or
                method == 'jane4' or
                method == 'riata-hgt'):
            output = parse_hgts(input_f=input_f,
                                method=method)
        elif method == 'consel':
            output = parse_consel(input_f=input_f)
        elif method == 'darkhorse':
            output = parse_darkhorse(input_f=input_f,
                                     output_fp=output_fp,
                                     low_lpi=low_lpi,
                                     high_lpi=high_lpi)
        elif method == 'hgtector':
            output = parse_hgtector(input_f=input_f)
        else:
            raise ValueError("Method is not supported: %s" % method)
        return output


@click.command()
@click.option('--hgt-results-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Output file containing HGT information')
@click.option('--method', required=True,
              type=click.Choice(['trex', 'ranger-dtl',
                                 'riata-hgt', 'consel',
                                 'darkhorse', 'wn-svm',
                                 'genemark', 'hgtector',
                                 'distance-method', 'jane4',
                                 'tree-puzzle']),
              help='The method used for HGT detection')
@click.option('--darkhorse-low-lpi', type=float, default=0.0,
              show_default=True, required=False, help='Lower bound LPI score')
@click.option('--darkhorse-high-lpi', type=float, default=0.6,
              show_default=True, required=False, help='Upper bound LPI score')
@click.option('--darkhorse-output-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output all best hit IDs from DarkHorse summary')
def main(hgt_results_fp,
         method,
         darkhorse_low_lpi,
         darkhose_high_lpi,
         darkhorse_output_fp=None):
    """ Parsing functions for various HGT detection tool outputs.
    """
    output = parse_output(hgt_results_fp=hgt_results_fp,
                          method=method,
                          low_lpi=darkhorse_low_lpi,
                          high_lpi=darkhorse_high_lpi,
                          output_fp=darkhorse_output_fp)
    sys.stdout.write(output)


if __name__ == "__main__":
    main()
