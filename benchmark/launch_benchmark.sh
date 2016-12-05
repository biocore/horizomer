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
# to launch_software.sh for executing HGT detection software. For fields
# that are empty, None should be passed.

# working dir
working_dir=$(readlink -m $1)
# scripts dir
scripts_dir=$(readlink -m $2)
# species tree in Newick format
species_tree_fp=$3
# species genome in GenBank format (for compositional methods)
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
diamond_tabular_query_fp=${10}
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
init_command="${15}"
# DIAMOND NR
diamond_nr=${16}
# Number threads
threads=${17}
# Bash config file path (if None, default ~/.bash_profile)
bash_config="${18}"
# Launch on qsub cluster environment (true or false, if None, defaults to true)
qsub_env=${19}
# DarkHorse LPI upper bound
lpi_upper=${20}
# DarkHorse LPI lower bound
lpi_lower=${21}
# Parse HGTs for DarkHorse
parse_hgts=${22}

if [ "${init_command}" == "None" ]
then
    init_command="sleep 1"
fi

## Step 1: 
##    Align with DIAMOND query vs. NR
if [ "${diamond_tabular_query_fp}" == "None" ]
then
    if [ "$verbose" == "true" ]
    then
        echo "Running DIAMOND .."
    fi
    ## Build database if doesn't exist
    if [ "${diamond_nr_fp}" == "None" ]
    then
        diamond_nr_fp=${working_dir}/diamond/$(basename ${database_fp%.*})
        diamond makedb --in ${database_fp} -d ${diamond_nr_fp} --threads $threads
    fi
    ## Run DIAMOND
    filename=$(basename "${query_species_coding_seqs_fp}")
    diamond_output=${working_dir}/diamond/$filename
    diamond blastp --db ${diamond_nr_fp} \
                   --query ${query_species_coding_seqs_fp} \
                   --evalue 1e-5 \
                   --max-target-seqs 500 \
                   --threads ${threads} \
                   --daa ${diamond_output}.daa \
                   --sensitive
    # convert output to tab delimited format
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
    diamond_tabular_query_fp=${diamond_output}.m8
    if [ "$verbose" == "true" ]
    then
        echo "Done"
    fi
fi

# load submit_job function
. $scripts_dir/utils.sh

## Step 2:
##    Run DarkHorse and choose candidate reference genomes
##    (all genomes in DarkHorse output)
## Run DarkHorse
parse_hgts="false"
cmd="${init_command}; \
      bash ${scripts_dir}/run_darkhorse.sh ${nr_fp} \
                                           ${diamond_db_nr} \
                                           ${diamond_tabular_query_fp} \
                                           ${darkhorse_config_fp} \
                                           ${darkhorse_install_dp} \
                                           ${query_species_coding_seqs_fp} \
                                           ${working_dir} \
                                           ${threads} \
                                           ${verbose} \
                                           ${lpi_upper} \
                                           ${lpi_lower} \
                                           ${parse_hgts}"
submit_job "${cmd}" darkhorse

## Step 3:
##    Run PhyloPhlAn on candidate species genomes
if [ "${species_tree_fp}" == "None" ]
then
    continue
fi

## TODO Step 4:
##    Run Phylomizer to detect orthologous genes and build gene trees.
##    Use reference genomes from Step 2
if [ "${gene_tree_dir}" == "None" ]
then
    continue
fi

## Step 5: Launch all software
bash ${scripts_dir}/launch_software.sh ${working_dir} \
                                       ${scripts_dir} \
                                       ${species_tree_fp} \
                                       ${species_genome_fp} \
                                       ${species_model_fp} \
                                       ${query_species_coding_seqs_fp} \
                                       ${ref_species_coding_seqs_fp} \
                                       ${gene_tree_dir} \
                                       ${gene_msa_dir} \
                                       ${diamond_tabular_query_fp} \
                                       ${phylonet_install_dir} \
                                       ${jane_install_dir} \
                                       ${trex_install_dir} \
                                       ${verbose} \
                                       "${init_command}" \
                                       ${threads} \
                                       "${bash_config}" \
                                       ${qsub_env}





