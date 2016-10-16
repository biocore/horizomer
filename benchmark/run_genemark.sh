#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run GeneMark software
species_genome_fp=$1
output_fp=$2
stdout=$3
stderr=$4
scripts_dir=$5
genemark_install_dir=$6
working_dir=$7

## run GeneMarkS training (generate typical and atypical gene models)
mkdir -p "${working_dir}/input"
python ${scripts_dir}/reformat_input_gi.py --genbank-fp ${species_genome_fp} --output-dir "${working_dir}/input" --method 'egid'
PWD=$(pwd)
cd ${working_dir}
cp ${genemark_install_dir}/.gm_key ./
${genemark_install_dir}/gmsn.pl --combine --gm --clean --name id input/id.fna 1>>$stdout 2>>$stderr
${genemark_install_dir}/gmhmmp -r -m id_hmm_combined.mod -o GeneMark_output.txt input/id.fna 1>>$stdout 2>>$stderr
python ${scripts_dir}/parse_output_genemark.py --genbank-fp input/id.gbk --genemark-output-fp GeneMark_output.txt >> $output_fp
rm -rf input
rm id*
cd $PWD
