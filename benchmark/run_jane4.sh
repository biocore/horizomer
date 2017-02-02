#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run JANE 4 software
set -eu
source $(dirname "$0")/utils.sh
args=(
    gene_tree_dir
    species_tree_fp
    input_file_nex
    output_file
    output_fp
    scripts_dir
    jane_install_dir
    stdout
    stderr
    verbose
)
get_args "$@"

TIMEFORMAT='%U %R'
total_user_time_jane="0.0"
total_wall_time_jane="0.0"
i=0
printf "#JANE4\n" >> $output_fp
touch ${output_file%.*}.total_results.txt

# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $output_fp

    # JANE4
    # input conditions: requires NEXUS input file;
    # supports one-to-many mapping in both directions (ex. multiple genes per species)
    python ${scripts_dir}/reformat_input.py --method 'jane4' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nex
    TIME="$( time ($jane_install_dir/jane-cli.sh $input_file_nex 1>${output_file} 2>>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp ${output_file} --method 'jane4' >> $output_fp
    printf "\n" >> $output_fp
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_jane=$(echo $total_user_time_jane + $user_time | bc)
    total_wall_time_jane=$(echo $total_wall_time_jane + $wall_time | bc)
    echo "#!#Gene $i" >> ${output_file%.*}.total_results.txt
    cat $output_file >> ${output_file%.*}.total_results.txt
    rm $output_file
    i=$((i+1))
done

echo "Total wall time Jane 4: $total_wall_time_jane" >> $output_fp
echo "Total user time Jane 4: $total_user_time_jane" >> $output_fp
