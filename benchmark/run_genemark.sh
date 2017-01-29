#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run GeneMark software
set -eu
args=(
  working_dir
  species_genome_fp
  output_fp
  scripts_dir
  genemark_install_dir
  stdout
  stderr
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
