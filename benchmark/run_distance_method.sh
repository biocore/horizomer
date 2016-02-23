#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run the Distance Method algorithm

distance_method_install_dir=$1
query_species_coding_seqs_fp=$2
target_proteomes_dir=$3
stdout=$3
stderr=$4
working_dir=$5
threads=$6

python $distance_method_install_dir/distance_method.py $query_species_coding_seqs_fp \ 
													   $target_proteomes_dir \ 
													   $working_dir \ 
													   --threads $threads \
													   --alignment-method diamond