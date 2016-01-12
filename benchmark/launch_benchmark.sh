#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# This script will preprocess files (determine orthologous genes, compute
# their trees and MSAs, run BLASTP on the query genome) and passes this data
# to launch_software.sh for executing HGT detection software

# usage: bash launch_benchmark.sh working_dir scripts_dir species_tree_fp \
#             species_genome_fp species_model_fp query_species_coding_seqs_fp \
#             ref_species_coding_seqs_fp gene_tree_dir gene_msa_dir \
#             phylonet_install_dir jane_install_dir trex_install_dir \
#             verbose_str

# working dir
working_dir=$(readlink -m $1)
# scripts dir
scripts_dir=$(readlink -m $2)
# species tree in Newick format
species_tree_fp=$3
# species raw genome in FASTA format
species_genome_fp=$4
# species HMM model (produced by GeneMarkS)
species_model_fp=$5
# query species protein coding sequences in FASTA format
query_species_coding_seqs_fp=$6
# reference species protein coding sequences in FASTA format
ref_species_coding_seqs_fp=$7
# gene trees in Newick format
gene_tree_dir=$8
# gene multiple sequence alignment dir
gene_msa_dir=$9
# tabular DIAMOND alignments of query genome
diamond_query_tabular=${10}
# PhyloNet install dir
phylonet_install_dir=${11}
# Jane 4 install dir
jane_install_dir=${12}
# T-REX install dir
trex_install_dir=${13}
# Verbose string 'true' or 'false'
verbose=${14}
# Initial command that precedes call to software
# (example choosing virtualenv to workon)
init_command=${15}
# DIAMOND NR
diamond_nr=${16}
# Number threads
threads=${17}

if [ "$verbose" == "true" ]
then
    echo "working dir: $working_dir"
    echo "scripts_dir: $scripts_dir"
    echo "species tree: $species_tree_fp"
    echo "species genome: $species_genome_fp"
    echo "HMM model: $species_model_fp"
    echo "query genomes: $query_species_coding_seqs_fp"
    echo "species proteome: $ref_species_coding_seqs_fp"
    echo "gene trees: $gene_tree_dir"
    echo "gene MSAs: $gene_msa_dir"
    echo "blast nr: ${blast_nr}"
    echo "PhyloNet install dir: $phylonet_install_dir"
    echo "Jane 4 install dir: $jane_install_dir"
    echo "T-REX install dir: $trex_install_dir"
fi

mkdir -p "${working_dir}/diamond"
filename=$(basename "${query_species_coding_seqs_fp}")
diamond_output=${working_dir}/diamond/$filename
# launch DIAMOND if alignments don't exist for genome
if [ ! -f "${diamond_output}"]
then
    diamond blastp --db ${diamond_nr} \
                   --query ${query_species_coding_seqs_fp} \
                   --evalue 1e-5 \
                   --max-target-seqs 500 \
                   --threads ${threads} \
                   --daa ${diamond_output}.daa \
                   --sensitive
    # convert output to tab delimited format
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
fi

# build HMM model if doesn't exist for genome


# launch all software
bash ${scripts_dir}/launch_software.sh ${working_dir} \
                                       ${scripts_dir} \
                                       ${species_tree_fp} \
                                       ${species_genome_fp} \
                                       ${species_model_fp} \
                                       ${query_species_coding_seqs_fp} \
                                       ${ref_species_coding_seqs_fp} \
                                       ${gene_tree_dir} \
                                       ${gene_msa_dir} \
                                       ${diamond_output}.m8 \
                                       ${phylonet_install_dir} \
                                       ${jane_install_dir} \
                                       ${trex_install_dir} \
                                       ${verbose} \
                                       ${init_command} \
                                       ${threads}





