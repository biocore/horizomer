#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run GeneMark software
set -eu
source $(dirname "$0")/utils.sh
args=(
    working_dir
    species_genome_fp
    output_fp
    scripts_dir
    genemark_install_dir
    stdout
    stderr
)
get_args "$@"

## run GeneMarkS training (generate typical and atypical gene models)
TIMEFORMAT='%U %R'
mkdir -p "${working_dir}/input.genemark"
mkdir -p "${working_dir}/output.genemark"
python ${scripts_dir}/reformat_input.py --genbank-fp ${species_genome_fp} --output-dir "${working_dir}/input.genemark" --method 'genemark'
PWD=$(pwd)
cd ${working_dir}/output.genemark
cp ${genemark_install_dir}/.gm_key ./
TIME="$( time (${genemark_install_dir}/gmsn.pl --verbose --combine --gm --clean --name id ../input.genemark/id.fna 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $TIME | awk '{print $1}')
wall_time=$(echo $TIME | awk '{print $2}')
TIME="$( time (${genemark_install_dir}/gmhmmp -v -r -m id_hmm_combined.mod -o GeneMark_output.txt ../input.genemark/id.fna 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $user_time + $(echo $TIME | awk '{print $1}') | bc | awk '$0+=$3')
wall_time=$(echo $wall_time + $(echo $TIME | awk '{print $2}') | bc | awk '$0+=$3')
python ${scripts_dir}/parse_output.py --genbank-fp ../input.genemark/id.gbk --hgt-results-fp GeneMark_output.txt --method 'genemark' >> $output_fp
# Clean up
rm .gm_key
cd $PWD

echo "Total user time GeneMark: ${user_time}" >> $stderr
echo "Total wall time GeneMark: ${wall_time}" >> $stderr
