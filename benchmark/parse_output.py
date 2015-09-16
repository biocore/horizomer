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


def parse_trex(input_f):
    """ Parse output of T-REX version 3.6.

    Parameters
    ----------
    input_f: string
        file descriptor for T-REX output results

    Returns
    -------
    number_of_hgts: string
        number of HGTs reported by a tool, or NaN if an error occurred
    """
    string = 'hgt : number of HGT(s) found = '
    for line in input_f:
        if string in line:
            return line.split(string)[1].strip()
    return "NaN"


def parse_rangerdtl(input_f):
    """ Parse output of RANGER-DTL version 1.0.

    Parameters
    ----------
    input_f: string
        file descriptor for RANGER-DTL output results

    Returns
    -------
    number_of_hgts: string
        number of HGTs reported by a tool, or NaN if an error occurred
    """
    string = "The minimum reconciliation cost is: "
    for line in input_f:
        if string in line:
            return line.split("Transfers: ")[1].split(",")[0]
    return "NaN"


def parse_riatahgt(input_f):
    """ Parse output of RIATA-HGT version 3.5.6.

    Parameters
    ----------
    input_f: string
        file descriptor for RIATA-HGT output results

    Returns
    -------
    number_of_hgts: string
        number of HGTs reported by a tool, or NaN if an error occurred
    """
    string = "There are "
    for line in input_f:
        if string in line:
            return line.split(string)[1].split(" component(s)")[0]
    return "NaN"


def parse_jane4(input_f):
    """ Parse output of Jane 4.

    Parameters
    ----------
    input_f: string
        file descriptor for RIATA-HGT output results

    Returns
    -------
    number_of_hgts: string
        number of HGTs reported by a tool, or NaN if an error occurred
    """
    string = "Host Switch: "
    for line in input_f:
        if string in line:
            return line.split(string)[1].strip()
    return "NaN"


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
    input_f.next()
    input_f.next()
    input_f.next()
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
        if method == 'ranger-dtl':
            output = parse_rangerdtl(input_f=input_f)
        elif method == 'trex':
            output = parse_trex(input_f=input_f)
        elif method == 'riata-hgt':
            output = parse_riatahgt(input_f=input_f)
        elif method == 'jane4':
            output = parse_jane4(input_f=input_f)
        elif method == 'consel':
            output = parse_consel(input_f=input_f)
            if output is None:
                output = "NaN"
            else:
                output = " ".join(output)

    sys.stdout.write(output)


if __name__ == "__main__":
    _main()
