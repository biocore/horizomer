#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run CONSEL software
gene_tree_dir=$1
output_fp=$2
verbose=$3
stdout=$4
stderr=$5
scripts_dir=$6
species_tree_fp=$7
input_file_nwk=$8
output_file=$9
gene_msa_dir=${10}
working_dir=${11}

TIMEFORMAT='%U %R'
total_user_time_consel="0.0"
total_wall_time_consel="0.0"
i=0
printf "y\n" > $working_dir/puzzle_cmd.txt
printf "#CONSEL\n" >> $output_fp

# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $output_fp

    # CONSEL (AU Test)
    # input conditions: matrix of the site-wise log-likelihoods
    # (a) if no MSA provided, CLUSTALW (align sequences)
    # (b) if MSA provided (ex. ALF), Fasta2Phylip.py
    # TREE-PUZZLE (reconstruct phylogenetic tree using maximum likelihood)
    # CONSEL (apply AU Test on matrix)
    gene_msa_fasta_fp=$gene_msa_dir/"MSA_${gene_number}_aa.fa"
    gene_msa_phylip_fp=$working_dir/"MSA_${gene_number}_aa.phy"
    python ${scripts_dir}/reformat_input.py --method 'tree-puzzle' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --gene-msa-fa-fp $gene_msa_fasta_fp \
                                            --output-tree-fp $input_file_nwk \
                                            --output-msa-phy-fp $gene_msa_phylip_fp
    puzzle -wsl $gene_msa_phylip_fp ${input_file_nwk} < $working_dir/puzzle_cmd.txt 1>>$stdout 2>>$stderr
    # makermt removes the .sitelh extension and writes to the edited file path
    # which would overwrite the Newick tree. Rename the input file to avoid this.
    mv ${input_file_nwk}.sitelh ${input_file_nwk}_puzzle.sitelh
    TIME="$( time (makermt --puzzle ${input_file_nwk}_puzzle.sitelh 1>>$stdout 2>>$stderr) 2>&1)"
    consel ${input_file_nwk}_puzzle 1>>$stdout 2>>$stderr
    catpv ${input_file_nwk}_puzzle.pv 1>$output_file 2>>$stderr
    python ${scripts_dir}/parse_output.py --hgt-results-fp $output_file --method 'consel' >> $output_fp
    printf "\n" >> $output_fp
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_consel=$(echo $total_user_time_consel + $user_time | bc)
    total_wall_time_consel=$(echo $total_wall_time_consel + $wall_time | bc)

    ## Clean up
    rm $output_file
    rm $gene_msa_phylip_fp
    rm ${input_file_nwk}_puzzle.pv
    rm ${input_file_nwk}_puzzle.sitelh
    i=$((i+1))
done

echo "Total wall time AU-Test: $total_wall_time_consel" >> $stderr
echo "Total user time AU-Test: $total_user_time_consel" >> $stderr