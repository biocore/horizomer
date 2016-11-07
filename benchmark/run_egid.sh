#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run EGID software
species_genome_fp=$1
output_fp=$2
stdout=$3
stderr=$4
scripts_dir=$5
egid_install_dir=$6
working_dir=$7

## run EGID
TIMEFORMAT='%U %R'
mkdir -p "${working_dir}/input.egid"
mkdir -p "${working_dir}/output.egid"
python ${scripts_dir}/reformat_input.py --genbank-fp ${species_genome_fp} --output-dir "${working_dir}/input.egid" --method 'egid'
PWD=$(pwd)
cd ${egid_install_dir}
TIME="$( time (./EGID "${working_dir}/input.egid" "${working_dir}/output.egid" 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $TIME | awk '{print $1}')                                                                                                                                                               
wall_time=$(echo $TIME | awk '{print $2}')
cd ${working_dir}
python ${scripts_dir}/parse_output.py --genbank-fp input.egid/id.gbk --hgt-results-fp output.egid/EGID_output.txt --method 'egid' >> $output_fp
# Clean up
rm input.egid/id.faa.mob
cd $PWD

echo "Total user time EGID: ${user_time}" >> $stderr
echo "Total wall time EGID: ${wall_time}" >> $stderr
