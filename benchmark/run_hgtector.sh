#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run HGTector software
args=(
  working_dir
  query_species_coding_seqs_fp
  hit_table_fp
  database_faa_fp
  database_dmnd_fp
  taxdump_dp
  gi_to_taxid_fp
  hgtector_config_fp
  threads
  taxid
  scripts_dir
  hgtector_install_dir
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

if [ "${gi_to_taxid_fp}" == "None" ]
then
    echo "GI to TaxID translation file is required."
    exit
fi

## Align with DIAMOND if alignments don't exist
if [ "${hit_table_fp}" == "None" ]
then
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
    # convert output to tab delimited format
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
    hit_table_fp=${diamond_output}.m8
fi

## Run HGTector
mkdir -p ${working_dir}/hgtector
mkdir -p ${working_dir}/hgtector/input
cp $query_species_coding_seqs_fp ${working_dir}/hgtector/input/id.faa
mkdir -p ${working_dir}/hgtector/presearch
ln -s $(readlink -f ${hit_table_fp}) \
   ${working_dir}/hgtector/presearch/id.m8

## create config.txt (if wasn't passed)
if [ "${hgtector_config_fp}" == "None" ]
then
    hgtector_config_fp=${working_dir}/hgtector/config.txt
    touch $hgtector_config_fp
    echo -e "title=trial_001\n\
    interactive=0\n\
    searchTool=DIAMOND\n\
    protdb=${database_dmnd_fp}.dmnd\n\
    taxdump=${taxdump_dp}\n\
    prot2taxid=${gi_to_taxid_fp}\n\
    preSearch=${working_dir}/hgtector/presearch\n\
    evalue=1e-20\n\
    identity=30\n\
    coverage=50\n\
    minSize=30\n\
    ignoreSubspecies=1\n\
    graphFp=1\n\
    exOutlier=2\n" > "$hgtector_config_fp"
    if [ "$taxid" != "None" ]
    then
        echo "selfTax=id:$taxid" >> "$hgtector_config_fp"
    fi
    if [ "$threads" != "None" ]
    then
        echo "threads=$threads" >> "$hgtector_config_fp"
    fi
else
    cp -f $hgtector_config_fp ${working_dir}/hgtector/config.txt
fi

perl ${hgtector_install_dir}/HGTector.pl ${working_dir}/hgtector

output_file=${working_dir}/hgtector/result/detail/id.txt
python ${scripts_dir}/parse_output.py --hgt-results-fp ${output_file} --method 'hgtector' >> $output_fp
