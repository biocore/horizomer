#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run DIAMOND to perform protein sequence similarity search
set -eu
args=(
  query_faa_fp
  database_faa_fp
  database_dmnd_fp
  output_hit_table
  working_dir
  scripts_dir
  verbose
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

if [ "$verbose" == "true" ]
then
    echo "Running DIAMOND .."
fi
mkdir -p ${working_dir}/diamond

# build DIAMOND database if doesn't exist
if [ "${database_dmnd_fp}" == "None" ]
then
    database_dmnd_fp=${working_dir}/diamond/$(basename ${database_faa_fp%.*})
    diamond makedb --in ${database_faa_fp} -d ${database_dmnd_fp} --threads $threads
fi

# run DIAMOND
filename=$(basename "${query_faa_fp}")
diamond_output=${working_dir}/diamond/$filename
diamond blastp --db ${database_dmnd_fp} \
               --query ${query_faa_fp} \
               --evalue 1e-5 \
               --max-target-seqs 500 \
               --threads ${threads} \
               --daa ${diamond_output}.daa \
               --sensitive

# convert output to tab delimited format
#   note: newer versions of DIAMOND can generate tabular output directly. The
#   current script is backward compatible with DIAMOND 7.10-.
diamond view --daa ${diamond_output}.daa -f tab -o ${output_hit_table}.m8
rm ${output_hit_table}.daa
hit_table_fp=${output_hit_table}.m8
if [ "$verbose" == "true" ]
then
    echo "Done"
fi
