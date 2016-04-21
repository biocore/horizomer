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
def _main(hgt_results_fp,
          method):
    """ Parsing functions for various HGT detection tool outputs.
    """
    with open(hgt_results_fp, 'U') as input_f:
        if (method == 'ranger-dtl' or
                method == 'trex' or
                method == 'jane4' or
                method == 'riata-hgt'):
            output = parse_hgts(input_f=input_f,
                                method=method)
        elif method == 'consel':
            output = parse_consel(input_f=input_f)
            if output is None:
                output = "NaN"
            else:
                output = " ".join(output)
        else:
            raise ValueError("Method is not supported: %s" % method)

    sys.stdout.write(output)


if __name__ == "__main__":
    _main()
