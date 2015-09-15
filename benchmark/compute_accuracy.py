# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Given known HGTs, compute precision, recall and F-score for various tools
=========================================================================
"""

import sys
import click
import re


def parse_expected_transfers(ground_truth_f):
    """ Parse ALF's log file.

    Parameters
    ----------
    ground_truth_f: string
        file descriptor to ALF's log file

    Returns
    -------
    expected_transfers: list of tuples
        list of transfers with each tuple representing (organism donor, gene
        donated, organism recipient, gene received)

    Notes
    -----
    Report horizontal gene transfers with information on organism donor,
    gene donated, organism recipient and gene received.
    """
    expected_transfers = []
    string = "lgt from organism "
    for line in ground_truth_f:
        if string in line:
            content = re.split(
                'lgt from organism | with gene | to organism |, now gene ',
                line.strip())
            expected_transfers.append(tuple(content[1:]))
    return expected_transfers


def parse_observed_transfers(observed_hgts_f,
                             pvalue_cutoff):
    """ Parse summary file of observed transfers for various tools.

    Parameters
    ----------
    observed_hgts_f: string
        file descriptor of observed transfers (output of launch_software.sh)
    pvalue_cutoff: float
        p-value threshold below which to consider two tree topologies
        incongruent (incongruence attributed to HGT)

    Returns
    -------
    observed_transfers: dict
        dictionary of tools' names (keys) and a list of horizontal gene
        transfers

    Notes
    -----
    Parse output of launch_software.sh. The outpur file is expected to
    follow the format:

    #number of HGTs detected
    #	gene ID	T-REX	RANGER-DTL	RIATA-HGT	Jane 4	Consel
    0	1000	1	1	1	1	0.96 0.04
    1	1001	0	0	0	0	0.00 0.00
    ..
    """
    tools = {}
    tools_id = []
    next(observed_hgts_f)
    for line in observed_hgts_f:
        line = line.strip().split('\t')
        if line[0].startswith('#'):
            for tool in line[2:]:
                tools_id.append(tool)
                if tool not in tools:
                    tools[tool] = []
            continue
        gene_id = line[1]
        for i, entry in enumerate(line[2:]):
            entry = entry.split()
            # consel output (2 p-values for AU Test)
            if len(entry) == 2:
                # species tree p-value
                if entry[0] != "0.00":
                    # gene tree p-value
                    pvalue = entry[1]
                    if pvalue != 'NaN':
                        if float(pvalue) <= pvalue_cutoff:
                            tools[tools_id[i]].append(gene_id)
            else:
                hgt_num = entry[0]
                if hgt_num != 'NaN':
                    if int(hgt_num) > 0:
                        tools[tools_id[i]].append(gene_id)
    return tools


def compute_accuracy(expected_transfers,
                     observed_transfers):
    """ Compute precision, recall and F-score for horizontally detected genes.

    Parameters
    ----------
    expected_transfers: list of tuples
        list of transfers with each tuple representing (organism donor, gene
        donated, organism recipient, gene received)
    observed_transfers: dict
        dictionary of tools' names (keys) and a list of horizontal gene
        transfers

    Returns
    -------
    tools_accuracy: dict
        dictionary of tuples for the number of true-positive (tp),
        false-positive (fp), false-negative (fn) reported HGTs and precision,
        recall and F-score
    """
    exp_s = set()
    obs_s = set()
    tools_accuracy = {}
    for tup in expected_transfers:
        exp_s.add(tup[1])
    sys.stdout.write("#expected HGTs: %s\n" % len(exp_s))
    sys.stdout.write("#tool\tTP\tFP\tFN\tprecision\trecall\tF-score\n")
    for tool in observed_transfers:
        obs_s = set(observed_transfers[tool])
        if not obs_s:
            continue
        tp = len(obs_s & exp_s)
        fp = len(obs_s - exp_s)
        fn = len(exp_s - obs_s)
        p = tp / float(tp + fp)
        r = tp / float(tp + fn)
        f = float(2 * p * r) / float(p + r)
        tools_accuracy[tool] = (tp, fp, fn, p, r, f)
    return tools_accuracy


@click.command()
@click.option('--ground-truth-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='logfile.txt from ALF simulations')
@click.option('--observed-hgts-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='output from launch_software.sh')
@click.option('--pvalue-cutoff', required=False, default=0.05,
              show_default=True, help='p-value threshold below which to '
                                      'consider two tree topologies '
                                      'incongruent (AU Test)')
def _main(ground_truth_fp,
          observed_hgts_fp,
          pvalue_cutoff):
    """ Compute precision, recall and F-score for observed HGT, loss and gain.
    """
    with open(ground_truth_fp, 'U') as ground_truth_f:
        expected_transfers = parse_expected_transfers(ground_truth_f)
    with open(observed_hgts_fp, 'U') as observed_hgts_f:
        observed_transfers = parse_observed_transfers(observed_hgts_f,
                                                      pvalue_cutoff)

    tools_accuracy = compute_accuracy(expected_transfers, observed_transfers)
    for result in tools_accuracy:
        tup = tools_accuracy[result]
        sys.stdout.write("%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n" % (
            result, tup[0], tup[1], tup[2], tup[3], tup[4], tup[5]))


if __name__ == "__main__":
    _main()
