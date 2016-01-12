#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run RANGER-DTL software
gene_tree_dir=$1
output_fp=$2
verbose=$3
stdout=$4
stderr=$5
scripts_dir=$6
species_tree_fp=$7
input_file_nwk=$8
output_file=$9

total_user_time_rangerdtl="0.0"
total_wall_time_rangerdtl="0.0"
i=0

printf "#RANGER\n" >> $output_fp
# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $output_fp

    python ${scripts_dir}/reformat_input.py --method 'ranger-dtl' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nwk
    TIME="$( time (ranger-dtl-U.linux -i $input_file_nwk -o ${output_file} 1>>$stdout 2>>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp ${output_file} --method 'ranger-dtl' >> $output_fp
    printf "\n" >> $output_fp
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_rangerdtl=$(echo $total_user_time_rangerdtl + $user_time | bc)
    total_wall_time_rangerdtl=$(echo $total_wall_time_rangerdtl + $wall_time | bc)
    rm ${output_file}
    i=$((i+1))
done

echo "Total time RANGER-DTL: $total_user_time_rangerdtl"