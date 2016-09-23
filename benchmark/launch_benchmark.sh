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
# species raw genome in FASTA format, None if does not exist
species_genome_fp=$4
# species HMM model (produced by GeneMarkS), None if does not exist
species_model_fp=$5
# query species protein coding sequences in FASTA format
query_species_coding_seqs_fp=$6
# gene trees in Newick format, None if does not exist
gene_tree_dir=$7
# gene multiple sequence alignment dir, None if does not exist
gene_msa_dir=$8
# tabular DIAMOND alignments of query genome, None if does not exist
diamond_tabular_query_fp=$9
# PhyloNet install dir
phylonet_install_dir=${10}
# Jane 4 install dir
jane_install_dir=${11}
# T-REX install dir
trex_install_dir=${12}
# Verbose string 'true' or 'false'
verbose=${13}
# Initial command that precedes call to software
# (example choosing virtualenv to workon)
init_command="${14}"
# Number threads
threads=${15}
# nr FASTA file
nr_fp=${16}
# directory that contains names.dmp and nodes.dmp
nr_tax_dp=${17}
# nr DIAMOND database, None if doesn't exist
diamond_nr_fp=${18}
# Bash config file path (if None, default ~/.bash_profile)
bash_config="${19}"
# Launch on qsub cluster environment (true or false, if None, defaults to true)
qsub_env=${20}
# Darkhorse config file
darkhorse_config_fp=${21}
# Darkhorse install directory
darkhorse_install_dp=${22}
# DarkHorse LPI upper bound, None if default
lpi_upper=${23}
# DarkHorse LPI lower bound, None if default
lpi_lower=${24}
# HGTector config file
hgtector_config_fp=${25}
# HGTector install directory
hgtector_install_dp=${26}
# GI-to-TaxID translation table
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip
gi_to_taxid_fp=${27}
# EGID install directory
egid_install_dp=${28}
# HGT summary file basename 
hgt_summary=$working_dir/"hgt_summary"
# qsub params
qsub="-q route -m abe -M jenya.kopylov@gmail.com -l nodes=1:ppn=${threads} -l walltime=230:00:00 -l pmem=10gb -l mem=20gb"

mkdir -p "${working_dir}"
if [ "${init_command}" == "None" ]
then
    init_command="sleep 1"
fi

if [ "${bash_config}" == "None" ]
then
    bash_config="~/.bash_profile"
fi
if [ "${qsub_env}" == "None" ]
then
    qsub_env="true"
fi

# load submit_job function
. ${scripts_dir}/utils.sh

if [ "$verbose" == "true" ]
then
    echo "INPUT PARAMETER SUMMARY"
    echo "-----------------------"
    echo "working_dir: ${working_dir}"
    echo "scripts_dir: ${scripts_dir}"
    echo "species_tree_fp: ${species_tree_fp}"
    echo "species_genome_fp: ${species_genome_fp}"
    echo "species_model_fp: ${species_model_fp}"
    echo "query_species_coding_seqs_fp: ${query_species_coding_seqs_fp}"
    echo "gene_tree_dir: ${gene_tree_dir}"
    echo "gene_msa_dir: ${gene_msa_dir}"
    echo "diamond_tabular_query_fp: ${diamond_tabular_query_fp}"
    echo "phylonet_install_dir: ${phylonet_install_dir}"
    echo "jane_install_dir: ${jane_install_dir}"
    echo "trex_install_dir: ${trex_install_dir}"
    echo "verbose: ${verbose}"
    echo "init_command: ${init_command}"
    echo "threads: $threads"
    echo "nr_fp: ${nr_fp}"
    echo "nr_tax_dp: ${nr_tax_dp}"
    echo "diamond_nr_fp: ${diamond_nr_fp}"
    echo "bash_config: ${bash_config}"
    echo "qsub_env: ${qsub_env}"
    echo "darkhorse_config_fp: ${darkhorse_config_fp}"
    echo "darkhorse_install_dp: ${darkhorse_install_dp}"
    echo "lpi_upper: ${lpi_upper}"
    echo "lpi_lower: ${lpi_lower}"
    echo "hgtector_config_fp: ${hgtector_config_fp}"
    echo "hgtector_install_dp: ${hgtector_install_dp}"
    echo "gi_to_taxid_fp: ${gi_to_taxid_fp}"
    echo "egid_install_dp: ${egid_install_dp}"
    echo "-----------------------"
    echo ""
fi

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
        echo "Running DIAMOND"
        echo "-----------------------"
    fi
    ## Build database if doesn't exist
    if [ "${diamond_nr_fp}" == "None" ]
    then
        if [ "$verbose" == "true" ]
        then
          echo "diamond makedb"
        fi
        mkdir ${working_dir}/diamond
        diamond_nr_fp=${working_dir}/diamond/$(basename ${nr_fp%.*})
        diamond makedb --in ${nr_fp} -d ${diamond_nr_fp} --threads $threads
    fi
    ## Run DIAMOND
    if [ "$verbose" == "true" ]
    then
      echo "diamond blastp"
    fi
    filename=$(basename "${query_species_coding_seqs_fp}")
    diamond_output=${working_dir}/diamond/$filename
    diamond blastp --db ${diamond_nr_fp} \
                   --query ${query_species_coding_seqs_fp} \
                   --evalue 1e-5 \
                   --max-target-seqs 500 \
                   --threads ${threads} \
                   --daa ${diamond_output}.daa \
                   --sensitive
    # Convert output to tab delimited format
    if [ "$verbose" == "true" ]
    then
      echo "diamond view"
    fi
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
    diamond_tabular_query_fp=${diamond_output}.m8
    if [ "$verbose" == "true" ]
    then
        echo "Done"
        echo "-----------------------"
    fi
fi

# Load submit_job function
. $scripts_dir/utils.sh

## Step 2:
##    Run DarkHorse and choose candidate reference genomes
##    (all genomes in DarkHorse output)
## Run DarkHorse
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
                                           ${scripts_dir} \
                                           ${hgt_summary}.darkhorse.txt \
                                           ${hgt_summary}.candidate_genomes.darkhorse.txt"

submit_job "${cmd}" darkhorse

## TODO Step 3:
##    Prune complete species tree to include only species from Step 2
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
#bash ${scripts_dir}/launch_software.sh ${working_dir} \
#                                       ${scripts_dir} \
#                                       ${species_tree_fp} \
#                                       ${species_genome_fp} \
#                                       ${species_model_fp} \
#                                       ${query_species_coding_seqs_fp} \
#                                       ${gene_tree_dir} \
#                                       ${gene_msa_dir} \
#                                       ${diamond_tabular_query_fp} \
#                                       ${phylonet_install_dir} \
#                                       ${jane_install_dir} \
#                                       ${trex_install_dir} \
#                                       ${verbose} \
#                                       "${init_command}" \
#                                       ${threads} \
#                                       ${nr_fp} \
#                                       ${nr_tax_dp} \
#                                       ${diamond_nr_fp} \
#                                       "${bash_config}" \
#                                       ${qsub_env} \
#                                       ${hgtector_config_fp} \
#                                       ${hgtector_install_dp} \
#                                       ${gi_to_taxid_fp} \
#                                       ${egid_install_dp} \
#                                       ${hgt_summary} \
#                                       ${qsub}





