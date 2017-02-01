#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DarkHorse software
set -eu
source $(dirname "$0")/utils.sh
args=(
    hit_table_fp
    darkhorse_config_fp
    darkhose_install_dir
    query_species_coding_seqs_fp
    working_dir
    verbose
    lpi_upper
    lpi_lower
    parse_hgts
    scripts_dir
    output_fp
)
get_args "$@"

if [ "$verbose" == "true" ]
then
    echo "Parameters:"
    for i in ${args[@]}
    do
      echo "$i = ${!i}"
    done
fi

# set default LPI upper and lower bounds
# see http://darkhorse.ucsd.edu/tutorial.html
if [ "${lpi_upper}" == "None" ]
then
    lpi_upper="0.6"
fi
if [ "${lpi_lower}" == "None" ]
then
    lpi_lower="0.2"
fi

# select list of species for "exclude list template"
# run DarkHorse
mkdir -p "${working_dir}/darkhorse"
cmd="perl ${darkhose_install_dir}/bin/darkhorse.pl -c ${darkhorse_config_fp} \
                                                   -t ${hit_table_fp} \
                                                   -g ${query_species_coding_seqs_fp} \
                                                   -e ${darkhose_install_dir}/templates/exclude_list_template"
if [ "$verbose" == "true" ]
then
    echo "Running DarkHorse .."
    echo "command: $cmd"
fi
# create and enter DarkHorse working directory
# (DarkHorse does not currently support writing to defined output directory)
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
