#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# purpose: launch HGT software on testing data sets, reformat input trees to
#          follow input format for each tool and parse output files to report
#          standardized statistics (number of HGTs, donors, recipients, gains
#          and losses)

# working dir
working_dir=$1
echo "working dir: $working_dir"
# scripts dir
scripts_dir=$2
echo "scripts_dir: $scripts_dir"
# species tree in Newick format
species_tree_fp=$3
echo "species tree: $species_tree_fp"
# species raw genome in FASTA format
species_genome_fp=$4
echo "species genome: $species_genome_fp"
# species HMM model (produced by GeneMarkS)
species_model_fp=$5
echo "HMM model: $species_model_fp"
# query species protein coding sequences in FASTA format
query_species_coding_seqs_fp=$6
echo "query genomes: $query_species_coding_seqs_fp"
# reference species protein coding sequences in FASTA format
ref_species_coding_seqs_fp=$7
echo "species proteome: $ref_species_coding_seqs_fp"
# gene trees in Newick format
gene_tree_dir=$8
echo "gene trees: $gene_tree_dir"
# gene multiple sequence alignment dir
gene_msa_dir=$9
echo "gene MSAs: $gene_msa_dir"
# PhyloNet install dir
phylonet_install_dir=${10}
echo "PhyloNet install dir: $phylonet_install_dir"
# Jane 4 install dir
jane_install_dir=${11}
echo "Jane 4 install dir: $jane_install_dir"
# T-REX install dir
trex_install_dir=${12}
echo "T-REX install dir: $trex_install_dir"


TIMEFORMAT='%U %R'
base_input_file_nwk="input_tree.nwk"
base_input_file_nex="input_tree.nex"
base_input_file_phy="input_msa.phy"
base_output_file="output_file.txt"
hgt_summary_file=$working_dir/"hgt_summary.txt"
input_file_nwk=$working_dir/$base_input_file_nwk
input_file_nex=$working_dir/$base_input_file_nex
output_file=$working_dir/$base_output_file
input_msa_phy=$working_dir/$base_input_file_phy
stderr=$working_dir/"stderr.txt"
stdout=$working_dir/"stdout.txt"

if [ ! -d "${working_dir}" ]; then
    mkdir $working_dir
fi

total_user_time_trex="0.0"
total_wall_time_trex="0.0"
total_user_time_rangerdtl="0.0"
total_wall_time_rangerdtl="0.0"
total_user_time_riatahgt="0.0"
total_wall_time_riatahgt="0.0"
total_user_time_jane="0.0"
total_wall_time_jane="0.0"
total_user_time_consel="0.0"
total_wall_time_consel="0.0"

printf "y\n" > $working_dir/puzzle_cmd.txt

i=0
printf "#number of HGTs detected\n" > $hgt_summary_file
printf "#\tgene ID\tT-REX\tRANGER-DTL\tRIATA-HGT\tJane 4\tConsel\n" >> $hgt_summary_file

# search for HGTs in each gene tree
for gene_tree in $gene_tree_dir/*.nwk
do
    gene_tree_file=$(basename $gene_tree)
    gene_number=$(echo $gene_tree_file | sed 's/[^0-9]*//g')
    printf "$i\t$gene_number\t" >> $hgt_summary_file
    echo "Processing gene $i .."
    # T-REX
    echo -e "\tT-REX"
    python ${scripts_dir}/reformat_input.py --method 'trex' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nwk
    cp $input_file_nwk $trex_install_dir
    TIME="$( time (cd $trex_install_dir; ./hgt3.4 -inputfile=$base_input_file_nwk 1>$stdout 2>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp $stdout --method 'trex' >> $hgt_summary_file
    printf "\t" >> $hgt_summary_file
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_trex=$(python -c "print $total_user_time_trex + $user_time")
    total_wall_time_trex=$(python -c "print $total_wall_time_trex + $wall_time")
    rm $stdout

    # RANGER-DTL
    echo -e "\tRANGER-DTL"
    python ${scripts_dir}/reformat_input.py --method 'ranger-dtl' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nwk
    TIME="$( time (ranger-dtl-U.linux -i $input_file_nwk -o $output_file 1>$stdout 2>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp $output_file --method 'ranger-dtl' >> $hgt_summary_file
    printf "\t" >> $hgt_summary_file
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_rangerdtl=$(python -c "print $total_user_time_rangerdtl + $user_time")
    total_wall_time_rangerdtl=$(python -c "print $total_wall_time_rangerdtl + $wall_time")
    rm $output_file

    # RIATA-HGT (in PhyloNet)
    echo -e "\tRIATA-HGT"
    python ${scripts_dir}/reformat_input.py --method 'riata-hgt' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nex
    TIME="$( time (java -jar $phylonet_install_dir/PhyloNet_3.5.6.jar $input_file_nex 1>$output_file 2>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp $output_file --method 'riata-hgt' >> $hgt_summary_file
    printf "\t" >> $hgt_summary_file
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_riatahgt=$(python -c "print $total_user_time_riatahgt + $user_time")
    total_wall_time_riatahgt=$(python -c "print $total_wall_time_riatahgt + $wall_time")
    rm $output_file

    # JANE4
    # input conditions: requires NEXUS input file;
    # supports one-to-many mapping in both directions (ex. multiple genes per
    # species)
    echo -e "\tJane 4"
    python ${scripts_dir}/reformat_input.py --method 'jane4' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --output-tree-fp $input_file_nex
    TIME="$( time ($jane_install_dir/jane-cli.sh $input_file_nex 1>$output_file 2>$stderr) 2>&1)"
    python ${scripts_dir}/parse_output.py --hgt-results-fp $output_file --method 'jane4' >> $hgt_summary_file
    printf "\t" >> $hgt_summary_file
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_jane=$(python -c "print $total_user_time_jane + $user_time")
    total_wall_time_jane=$(python -c "print $total_wall_time_jane + $wall_time")
    rm $output_file

    # CONSEL (AU Test)
    # input conditions: matrix of the site-wise log-likelihoods
    # (a) if no MSA provided, CLUSTALW (align sequences)
    # (b) if MSA provided (ex. ALF), Fasta2Phylip.py
    # TREE-PUZZLE (reconstruct phylogenetic tree using maximum likelihood)
    # CONSEL (apply AU Test on matrix)
    echo -e "\tTree-Puzzle and CONSEL"
    gene_msa_fasta_fp=$gene_msa_dir/"MSA_${gene_number}_aa.fa"
    gene_msa_phylip_fp=$working_dir/"MSA_${gene_number}_aa.phy"
    python ${scripts_dir}/reformat_input.py --method 'tree-puzzle' \
                                            --gene-tree-fp $gene_tree \
                                            --species-tree-fp $species_tree_fp \
                                            --gene-msa-fa-fp $gene_msa_fasta_fp \
                                            --output-tree-fp $input_file_nwk \
                                            --output-msa-phy-fp $gene_msa_phylip_fp
    puzzle -wsl $gene_msa_phylip_fp $input_file_nwk < $working_dir/puzzle_cmd.txt 1>$stdout 2>$stderr
    # makermt removes the .sitelh extension and writes to the edited file path
    # which would overwrite the Newick tree. Rename the input file to avoid this.
    mv ${input_file_nwk}.sitelh ${input_file_nwk}_puzzle.sitelh
    TIME="$( time (makermt --puzzle ${input_file_nwk}_puzzle.sitelh 1>$stdout 2>$stderr) 2>&1)"
    consel ${input_file_nwk}_puzzle 1>$stdout 2>$stderr
    catpv ${input_file_nwk}_puzzle.pv 1>$output_file 2>$stderr
    python ${scripts_dir}/parse_output.py --hgt-results-fp $output_file --method 'consel' >> $hgt_summary_file
    printf "\n" >> $hgt_summary_file
    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time_consel=$(python -c "print $total_user_time_consel + $user_time")
    total_wall_time_consel=$(python -c "print $total_wall_time_consel + $wall_time")

    ## Clean up
    rm $gene_msa_phylip_fp
    rm ${input_file_nwk}_puzzle.pv
    rm ${input_file_nwk}_puzzle.sitelh
    i=$((i+1))
done

echo "Total time T-REX: $total_user_time_trex"
echo "Total time RANGER-DTL: $total_user_time_rangerdtl"
echo "Total time RIATA-HGT: $total_user_time_riatahgt"
echo "Total time Jane 4: $total_user_time_jane"
echo "Total time AU-Test: $total_user_time_consel"
