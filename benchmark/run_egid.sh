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
scripts_dir=$2
working_dir=$3
egid_install_dir=$4
output_fp=$5

## run EGID
mkdir -p "${working_dir}/input"
mkdir -p "${working_dir}/output"
python ${scripts_dir}/reformat_input_gi.py --genbank-fp ${species_genome_fp} --output-dir "${working_dir}/input" --method 'egid'
PWD=$(pwd)
cd ${egid_install_dir}
./EGID "${working_dir}/input" "${working_dir}/output"
cd $PWD
rm "${working_dir}/input/id.faa.mob"
python ${scripts_dir}/parse_output_gi.py --genbank-fp "${working_dir}/input/id.gbk" --gi-fp "${working_dir}/output/EGID_output.txt" >> $output_fp
