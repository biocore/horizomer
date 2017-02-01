#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run T-REX software
set -eu
source $(dirname "$0")/utils.sh
args=(
    gene_tree_dir
    species_tree_fp
    input_file_nwk
    output_fp
    scripts_dir
    trex_install_dir
    stdout
    stderr
    verbose
)
get_args "$@"

TIMEFORMAT='%U %R'
total_user_time_trex="0.0"
total_wall_time_trex="0.0"
i=0

printf "#TREX\n" >> $output_fp
touch ${stdout%.*}.total_results.txt
# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $output_fp

    python ${scripts_dir}/reformat_input.py --method 'trex' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nwk
    cp $input_file_nwk $trex_install_dir
    base_input_file_nwk=$(basename $input_file_nwk)
    TIME="$( time (cd $trex_install_dir; ./hgt3.4 -inputfile=${base_input_file_nwk} 1>$stdout 2>>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp $stdout --method 'trex' >> $output_fp
    echo "#!#Gene $i" >> ${stdout%.*}.total_results.txt
    cat $stdout >> ${stdout%.*}.total_results.txt
    printf "\n" >> $output_fp
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_trex=$(echo $total_user_time_trex + $user_time | bc)
    total_wall_time_trex=$(echo $total_wall_time_trex + $wall_time | bc)
    rm $stdout
    i=$((i+1))
done

echo "Total wall time T-REX: $total_wall_time_trex" >> $output_fp
echo "Total user time T-REX: $total_user_time_trex" >> $output_fp
