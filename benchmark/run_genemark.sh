#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run GeneMark software
species_model_fp=$1
output_file=$2
species_genome_fp=$3
stdout=$4
stderr=$5
working_dir=$6

## run GeneMarkS training (generate typical and atypical gene models)
mkdir -p "${working_dir}/hmm-models"
filename=$(basename "${species_genome_fp}")
hmm_output=${working_dir}/hmm-models/$filename
if [ "${species_model_fp}" == "None" ]
then
    gmsn.pl --combine --gm --clean --name ${species_model_fp} "${hmm_output%.*}"
    species_model_fp="${hmm_output%.*}_hmm_combined.mod"
fi

TIMEFORMAT='%U %R'
TIME="$( time (gmhmmp -r -m ${species_model_fp} -o $output_file ${species_genome_fp} 1>$stdout 2>>$stderr) 2>&1)"                                                                                              
user_time=$(echo $TIME | awk '{print $1;}')                                                                                                                                                               
wall_time=$(echo $TIME | awk '{print $2;}')

echo "Total user time GeneMark: ${user_time}" >> $stderr
echo "Total wall time GeneMark: ${wall_time}" >> $stderr