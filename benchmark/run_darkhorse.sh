#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DarkHorse software
nr_fp=$1
diamond_nr_fp=$2
diamond_tabular_query_fp=$3
darkhorse_config_fp=$4
darkhose_install_dir=$5
query_species_coding_seqs_fp=$6
working_dir=$7
threads=$8
verbose=$9
lpi_upper=${10}
lpi_lower=${11}
scripts_dir=${12}
hgt_fp=${13}
output_fp=${14}

printf "#DarkHorse\n" >> $output_fp

mkdir -p "${working_dir}/diamond"
if [ "$verbose" == "true" ]
then
    echo "Parameters:"
    echo "nr_fp = ${nr_fp}"
    echo "diamond_nr_fp = ${diamond_nr_fp}"
    echo "diamond_tabular_query_fp = ${diamond_tabular_query_fp}"
    echo "darkhorse_config_fp = ${darkhorse_config_fp}"
    echo "darkhose_install_dir = ${darkhose_install_dir}"
    echo "query_species_coding_seqs_fp = ${query_species_coding_seqs_fp}"
    echo "working_dir = ${working_dir}"
    echo "threads = ${threads}"
    echo "verbose = ${verbose}"
    echo "lpi_upper = ${lpi_upper}"
    echo "lpi_lower = ${lpi_lower}"
    echo "scripts_dir = ${scripts_dir}"
    echo "hgt_fp = ${hgt_fp}"
    echo "output_fp = ${output_fp}"
fi

## Align with DIAMOND if alignments don't exist
if [ "${diamond_tabular_query_fp}" == "None" ]
then
    if [ "$verbose" == "true" ]
    then
        echo "Running DIAMOND .."
    fi
    ## Build database if doesn't exist
    if [ "${diamond_nr_fp}" == "None" ]
    then
        diamond_nr_fp=${working_dir}/diamond/$(basename ${nr_fp%.*})
        diamond makedb --in ${nr_fp} -d ${diamond_nr_fp} --threads $threads
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
    # Convert output to tab delimited format
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
    diamond_tabular_query_fp=${diamond_output}.m8
    if [ "$verbose" == "true" ]
    then
        echo "Done"
    fi
fi

## Select list of species for "exclude list template"
## Run DarkHorse
if [ "$verbose" == "true" ]
then
    echo "Running DarkHorse .."
    echo "command: ${darkhose_install_dir}/bin/darkhorse.pl \
                                          -c ${darkhorse_config_fp} \
                                          -t ${diamond_tabular_query_fp} \
                                          -g ${query_species_coding_seqs_fp} \
                                          -e ${darkhose_install_dir}/templates/exclude_list_template"
fi
## Create DarkHorse working directory and cd into it (DarkHorse does not currently support
## writing to defined output directory)
mkdir -p "${working_dir}/darkhorse"
PWD=$(pwd)
cd ${working_dir}/darkhorse
perl ${darkhose_install_dir}/bin/darkhorse.pl \
                                          -c ${darkhorse_config_fp} \
                                          -t ${diamond_tabular_query_fp} \
                                          -g ${query_species_coding_seqs_fp} \
                                          -e ${darkhose_install_dir}/templates/exclude_list_template
cd $PWD
if [ "$verbose" == "true" ]
then
    echo "Done"
fi

## Parse HGTs and output tabular format with Gene Query ID, Best Hit ID,
## Best Hit Taxid, Best Hit Species, Best Hit Lineage, % id, % coverage
## and LPI score (must be between lpi_lower and lpi_upper)
if [ "$parse_hgts" == "true" ]
then

    python ${scripts_dir}/parse_output.py --hgt-results-fp ${working_dir}/darkhorse/calcs_*/*_smry \
                                          --method 'darkhorse' \
                                          --darkhorse-low-lpi ${lpi_lower} \
                                          --darkhorse-high-lpi ${lpi_upper} \
                                          --output-fp $output_fp >> ${hgt_fp}
fi

