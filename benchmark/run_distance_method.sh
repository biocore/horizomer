#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run the Distance Method algorithm
set -eu
source $(dirname "$0")/utils.sh
args=(
    working_dir
    query_species_coding_seqs_fp
    target_proteomes_dir
    threads
    distance_method_install_dir
    stdout
    stderr
)
get_args "$@"

python $distance_method_install_dir/distance_method.py $query_species_coding_seqs_fp \
                                                       $target_proteomes_dir \
                                                       $working_dir \
                                                       --threads $threads 1>$stdout 2>$stderr
