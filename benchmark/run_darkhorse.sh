#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DarkHorse software
args=(
  database_faa_fp
  database_dmnd_fp
  diamond_tabular_query_fp
  darkhorse_config_fp
  darkhose_install_dir
  query_species_coding_seqs_fp
  working_dir
  threads
  verbose
  lpi_upper
  lpi_lower
  parse_hgts
  scripts_dir
  output_fp
)
arg_str=$(IFS=,; echo "${args[*]/%/:}" | tr '_' '-')
TEMP=`getopt -o "" -l $arg_str -n "$0" -- "$@"`
eval set -- "$TEMP"
while true ; do
  case "$1" in
    --?*) eval $(echo ${1:2} | tr '-' '_')=$2 ; shift 2 ;;
    --) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
  esac
done

mkdir -p "${working_dir}/diamond"
if [ "$verbose" == "true" ]
then
    echo "Parameters:"
    for i in ${args[@]}
    do
      echo "$i = ${!i}"
    done
fi

## Default LPI upper bound (used by DarkHorse)
## See http://darkhorse.ucsd.edu/tutorial.html
if [ "${lpi_upper}" == "None" ]
then
    lpi_upper="0.6"
fi
if [ "${lpi_lower}" == "None" ]
then
    lpi_lower="0.2"
fi

## Align with DIAMOND if alignments don't exist
if [ "${diamond_tabular_query_fp}" == "None" ]
then
    if [ "$verbose" == "true" ]
    then
        echo "Running DIAMOND .."
    fi
    ## Build database if doesn't exist
    if [ "${database_dmnd_fp}" == "None" ]
    then
        database_dmnd_fp=${working_dir}/diamond/$(basename ${database_faa_fp%.*})
        diamond makedb --in ${database_faa_fp} -d ${database_dmnd_fp} --threads $threads
    fi
    ## Run DIAMOND
    filename=$(basename "${query_species_coding_seqs_fp}")
    diamond_output=${working_dir}/diamond/$filename
    diamond blastp --db ${database_dmnd_fp} \
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
cmd="perl ${darkhose_install_dir}/bin/darkhorse.pl -c ${darkhorse_config_fp} \
                                                   -t ${diamond_tabular_query_fp} \
                                                   -g ${query_species_coding_seqs_fp} \
                                                   -e ${darkhose_install_dir}/templates/exclude_list_template"
if [ "$verbose" == "true" ]
then
    echo "Running DarkHorse .."
    echo "command: $cmd"
fi
## Create DarkHorse working directory and cd into it (DarkHorse does not currently support
## writing to defined output directory)
mkdir -p "${working_dir}/darkhorse"
PWD=$(pwd)
cd ${working_dir}/darkhorse
$cmd
cd $PWD
if [ "$verbose" == "true" ]
then
    echo "Done"
fi

## TODO: Parse HGTs
if [ "$parse_hgts" == "true" ]
then
    python ${scripts_dir}/parse_output.py --hgt-results-fp ${working_dir}/darkhorse/calcs_*/*_smry \
                                          --method 'darkhorse' >> $output_fp
fi
