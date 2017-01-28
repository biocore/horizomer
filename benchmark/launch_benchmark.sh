#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# This script will preprocess files (determine orthologous genes, compute
# their trees and MSAs, run BLASTP on the query genome) and passes this data
# to launch_software.sh for executing HGT detection software. For fields
# that are empty, None should be passed.

# -e: script will exit if any command fails
# -u: force initialization of all variables
set -eu

# declare command-line arguments
args=(
  # working dir
  working_dir
  # scripts dir
  scripts_dir
  # species tree in Newick format
  species_tree_fp
  # species genome in GenBank format (for compositional methods)
  species_genome_fp
  # species HMM model (produced by GeneMarkS)
  species_model_fp
  # query species protein coding sequences in FASTA format
  query_species_coding_seqs_fp
  # reference species protein coding sequences in FASTA format
  ref_species_coding_seqs_fp
  # gene trees in Newick format
  gene_tree_dir
  # gene multiple sequence alignment dir
  gene_msa_dir
  # tabular DIAMOND alignments of query genome
  diamond_tabular_query_fp
  # reference protein sequence database compiled by DIAMOND
  database_dmnd_fp
  # reference protein sequence database in Fasta format
  database_faa_fp
  # PhyloNet install dir
  phylonet_install_dir
  # Jane 4 install dir
  jane_install_dir
  # T-REX install dir
  trex_install_dir
  # Verbose string 'true' or 'false'
  verbose
  # Initial command that precedes call to software
  # (example choosing virtualenv to workon)
  init_command
  # Number threads
  threads
  # Bash config file path (if None, default ~/.bash_profile)
  bash_config
  # Launch on qsub cluster environment (true or false, if None, defaults to true)
  qsub_env
  # DarkHorse LPI upper bound
  lpi_upper
  # DarkHorse LPI lower bound
  lpi_lower
  # Parse HGTs for DarkHorse
  parse_hgts
)

# convert arguments to --long-options
arg_str=$(IFS=,; echo "${args[*]/%/:}" | tr '_' '-')

# use GNU getopt to retrieve arguments
TEMP=`getopt -o "" -l $arg_str -n "$0" -- "$@"`
eval set -- "$TEMP"
while true ; do
  case "$1" in
    --?*) eval $(echo ${1:2} | tr '-' '_')=$2 ; shift 2 ;;
    --) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
  esac
done

# manipulate arguments
working_dir=$(readlink -m $working_dir)
mkdir -p $working_dir
scripts_dir=$(readlink -m $scripts_dir)
init_command="$init_command"
if [ "${init_command}" == "None" ]
then
    init_command="sleep 1"
fi
bash_config="$bash_config"

## Step 1:
##    Align with DIAMOND query vs. reference database (e.g., NCBI nr)
if [ "${diamond_tabular_query_fp}" == "None" ]
then
    if [ "$verbose" == "true" ]
    then
        echo "Running DIAMOND .."
    fi
    mkdir -p ${working_dir}/diamond
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
    # Note: newer versions of DIAMOND can generate tabular output directly.
    diamond view --daa ${diamond_output}.daa -f tab -o ${diamond_output}.m8
    diamond_tabular_query_fp=${diamond_output}.m8
    if [ "$verbose" == "true" ]
    then
        echo "Done"
    fi
fi

# load submit_job function
. $scripts_dir/utils.sh

## Step 2:
##    Run DarkHorse and choose candidate reference genomes
##    (all genomes in DarkHorse output)
## Run DarkHorse
parse_hgts="false"
cmd="${init_command}; \
      bash ${scripts_dir}/run_darkhorse.sh ${nr_fp} \
                                           ${database_dmnd_fp} \
                                           ${diamond_tabular_query_fp} \
                                           ${darkhorse_config_fp} \
                                           ${darkhorse_install_dp} \
                                           ${query_species_coding_seqs_fp} \
                                           ${working_dir} \
                                           ${threads} \
                                           ${verbose} \
                                           ${lpi_upper} \
                                           ${lpi_lower} \
                                           ${parse_hgts}"
submit_job "${cmd}" darkhorse

## Step 3:
##    Run PhyloPhlAn on candidate species genomes
if [ "${species_tree_fp}" == "None" ]
then
    continue
fi

## TODO Step 4:
##    Run Phylomizer to detect orthologous genes and build gene trees.
##    Use reference genomes from Step 2
if [ "${gene_tree_dir}" == "None" ]
then
    continue
fi

## Step 5: Launch all software
bash ${scripts_dir}/launch_software.sh ${working_dir} \
                                       ${scripts_dir} \
                                       ${species_tree_fp} \
                                       ${species_genome_fp} \
                                       ${species_model_fp} \
                                       ${query_species_coding_seqs_fp} \
                                       ${ref_species_coding_seqs_fp} \
                                       ${gene_tree_dir} \
                                       ${gene_msa_dir} \
                                       ${diamond_tabular_query_fp} \
                                       ${phylonet_install_dir} \
                                       ${jane_install_dir} \
                                       ${trex_install_dir} \
                                       ${verbose} \
                                       "${init_command}" \
                                       ${threads} \
                                       "${bash_config}" \
                                       ${qsub_env}
