#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run RIATA-HGT software
gene_tree_dir=$1
output_fp=$2
verbose=$3
stdout=$4
stderr=$5
scripts_dir=$6
species_tree_fp=$7
input_file_nex=$8
output_file=$9
phylonet_install_dir=${10}

TIMEFORMAT='%U %R'
total_user_time_riatahgt="0.0"
total_wall_time_riatahgt="0.0"
i=0
printf "#RIATAHGT\n" >> $output_fp
touch ${output_file%.*}.total_results.txt

# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $output_fp

    python ${scripts_dir}/reformat_input.py --method 'riata-hgt' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nex
    TIME="$( time (java -jar $phylonet_install_dir/PhyloNet_3.5.7.jar $input_file_nex 1>$output_file 2>>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp ${output_file} --method 'riata-hgt' >> $output_fp
    printf "\n" >> $output_fp
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_riatahgt=$(echo $total_user_time_riatahgt + $user_time | bc)
    total_wall_time_riatahgt=$(echo $total_user_time_riatahgt + $user_time | bc)
    echo "#!#Gene $i" >> ${output_file%.*}.total_results.txt
    cat $output_file >> ${output_file%.*}.total_results.txt
    rm $output_file
    i=$((i+1))
done

echo "Total wall time RIATA-HGT: $total_wall_time_riatahgt" >> $output_fp
echo "Total user time RIATA-HGT: $total_user_time_riatahgt" >> $output_fp