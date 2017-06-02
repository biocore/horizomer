#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DIAMOND to perform protein sequence similarity search
# set -eu
source $(dirname "$0")/utils.sh
args=(
    query_faa_fp
    database_faa_fp
    database_dmnd_fp
    evalue
    pident
    qcovs
    max_target_seqs
    output_hit_table
    working_dir
    scripts_dir
    threads
    verbose
)
get_args "$@"

$verbose && echo "Running DIAMOND .."
mkdir -p ${working_dir}/diamond

# build DIAMOND database if doesn't exist
if [ "${database_dmnd_fp}" == "None" ]
then
    database_dmnd_fp=${working_dir}/diamond/$(basename ${database_faa_fp%.*})
    diamond makedb --in ${database_faa_fp} \
                   --db ${database_dmnd_fp} \
                   --threads $threads
fi

# run DIAMOND
filename=$(basename "${query_faa_fp}")
diamond_output=${working_dir}/diamond/$filename
diamond blastp --db ${database_dmnd_fp} \
               --query ${query_faa_fp} \
               --evalue ${evalue} \
               --id ${pident} \
               --query-cover ${qcovs} \
               --max-target-seqs ${max_target_seqs} \
               --threads ${threads} \
               --out ${output_hit_table} \
               --sensitive

$verbose && echo "Done"
