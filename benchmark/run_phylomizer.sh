#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run Phylomizer to build phylogenetic trees of individual gene families
# notes: the official "KMM" protocol is used

set -eu
source $(dirname "$0")/utils.sh
args=(
    gene_fa_dir  # protein sequences in multi-fasta format (*.fa), one file per gene family
    gene_tree_dir  # output gene trees will be saved to this directory
    output_fp
    scripts_dir
    phylomizer_install_dir
    phylomizer_config_fp
    threads
    stdout
    stderr
    verbose
)
get_args "$@"

$verbose && echo "Running Phylomizer .."

# set up
PWD=$(pwd)
mkdir -p ${working_dir}/phylomizer
if [ "${phylomizer_config_fp}" == "None" ]
then
    phylomizer_config_fp=${working_dir}/phylomizer/config.txt
    # modified KMM protocol
    # minimum number of sequences: 10
    # only maximum likelihood, but not neighbor joining, is to be conducted
    echo "
    verbose             parameter    1
    residue_datatype    parameter    protein
    force_seed_sequence parameter    True
    alignment           mode         kalign muscle mafft
    consensus           mode         m_coffee
    trimming            mode         trimal
    both_direction      parameter    True
    min_seqs            parameter    10
    in_letter           parameter    U:B
    in_letter           parameter    O:Z
    muscle              binary
    muscle_params       parameter
    mafft               binary
    mafft_params        parameter    --auto
    kalign              binary
    kalign_params       parameter    -f fasta
    m_coffee            binary       t_coffee
    m_coffee_params     parameter    -n_core 1 -output fasta -quiet
    trimal              binary
    trimal_params       parameter    -phylip -gt 0.1
    trimal_compare      parameter    -ct 0.1667
    readal              binary
    tree                mode         phyml
    evol_models         parameter    JTT WAG MtREV VT LG Blosum62 Dayhoff
    numb_models         parameter    2
    tree_approach       mode         ml
    ml                  parameter    -b -2 -o tlr
    phyml               binary
    phyml_params        parameter    -d aa -f e -v e -a e -c 4 --no_memory_check
    " > ${phylomizer_config_fp}
fi
$verbose && echo "Configuration file:"$'\n'"  ${phylomizer_config_fp}"

# command
cmd="python ${phylomizer_install_dir}/source/phylomizer.py -i input.fa --steps alignments trees -c ${phylomizer_config_fp} -o ."
$verbose && echo "Command:"$'\n'"  $cmd"

TIMEFORMAT='%U %R'
total_user_time="0.0"
total_wall_time="0.0"

printf "#Phylomizer\n" >> $output_fp

i=0
cd ${working_dir}/phylomizer
for gene_fa in $gene_fa_dir/*.fa
do
    # set up
    gene=$(basename $gene_fa)
    gene=${gene%.fa}
    mkdir -p $gene
    cd $gene
    ln -s ${gene_fa} input.fa

    # run Phylomizer and record time
    TIME="$( time ($cmd 1>$stdout 2>>$stderr) 2>&1)"

    # report the best aa substitution model and the -log likelihood of the tree
    best_tree=$(head -n1 input.tree.phyml.rank.ml)
    best_model=$(echo ${best_tree} | cut -f1)
    printf "$i\t$gene\t${best_tree}\n" >> $output_fp

    # copy the best tree to the gene tree directory
    cp input.tree.phyml.ml.${best_model}.nw ${gene_tree_dir}/$gene.nwk

    user_time=$(echo $TIME | awk '{print $1;}')
    wall_time=$(echo $TIME | awk '{print $2;}')
    total_user_time=$(echo $total_user_time + $user_time | bc)
    total_wall_time=$(echo $total_wall_time + $wall_time | bc)

    cd ..
    i=$((i+1))
done

echo "Total wall time Phylomizer: $total_wall_time" >> $output_fp
echo "Total user time Phylomizer: $total_user_time" >> $output_fp

cd $PWD

$verbose && echo "Done"
