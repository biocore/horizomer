#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DarkHorse software
database_fp=$1
diamond_nr_fp=$2
diamond_tabular_query_fp=$3
darkhorse_config_fp=$4
darkhose_install_dir=$5
query_species_coding_seqs_fp=$6
working_dir=$7
threads=$8

mkdir -p "${working_dir}/diamond"

## Align with DIAMOND if alignments don't exist
if [ "${diamond_tabular_query_fp}" == "None" ]
then
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
fi

## Select list of species for "exclude list template"
## Run DarkHorse
perl ${darkhose_install_dir}/darkhorse.pl -c ${darkhorse_config_fp} \ 
                                          -t ${diamond_tabular_query_fp} \ 
                                          -g ${query_species_coding_seqs_fp} \ 
                                          -e ${darkhose_install_dir}/templates/exclude_list_template


