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

# load utilities
source $(dirname "$0")/utils.sh

# declare command-line arguments
args=(
    ## program runtime behavior
    # working dir
    working_dir
    # scripts dir
    scripts_dir
    # verbose screen output (true or false)
    verbose
    # initial command that precedes call to software
    # (e.g., choosing virtualenv to work on)
    init_command
    # number of threads
    threads
    # Bash config file path (if None, default ~/.bash_profile)
    bash_config
    # launch on qsub cluster environment (true or false, if None, defaults to true)
    qsub_env

    ## input files and directories
    # species tree in Newick format
    species_tree_fp
    # species genome in GenBank format
    species_genbank_fp
    # species HMM model (produced by GeneMarkS)
    species_model_fp
    # query species protein coding sequences in FASTA format
    species_faa_fp
    # reference species protein coding sequences in FASTA format
    ref_species_faa_fp
    # gene trees in Newick format
    gene_tree_dir
    # gene multiple sequence alignment dir
    gene_msa_dir
    # sequence similarity search hit table in standard tabular format
    # (e.g., BLAST -outfmt 6 or DIAMOND tab)
    hit_table_fp

    ## reference databases
    # reference protein sequence database compiled by DIAMOND
    database_dmnd_fp
    # reference protein sequence database in Fasta format
    database_faa_fp

    ## application installation directories
    # Darkhorse
    darkhorse_install_dir
    # PhyloNet
    phylonet_install_dir
    # Jane 4
    jane_install_dir
    # T-REX
    trex_install_dir

    ## parameters for specific applications
    # DarkHorse LPI upper bound
    lpi_upper
    # DarkHorse LPI lower bound
    lpi_lower
    # parse HGTs for DarkHorse
    parse_hgts
)
get_args "$@"

# manipulate arguments
mkdir -p $working_dir
init_command="$init_command"
if [ "${init_command}" == "None" ]
then
    init_command="sleep 1"
fi
bash_config="$bash_config"

## Step 1:
##    Align with DIAMOND query vs. reference database (e.g., NCBI nr)
if [ "${hit_table_fp}" == "None" ]
then
    filename=$(basename "${species_faa_fp}")
    diamond_output=${working_dir}/diamond/$filename
    bash ${scripts_dir}/run_diamond.sh --query-faa-fp ${species_faa_fp} \
                                       --database-faa-fp ${database_faa_fp} \
                                       --database-dmnd-fp ${database_dmnd_fp} \
                                       --output-hit-table ${diamond_output} \
                                       --working-dir ${working_dir} \
                                       --scripts-dir ${scripts_dir} \
                                       --threads ${threads} \
                                       --verbose ${verbose}
    hit_table_fp=${diamond_output}.m8
fi

# load submit_job function
. $scripts_dir/utils.sh

if [ "${init_command}" != "None" ]
then
    ${init_command}
fi

## Step 2:
##    Run DarkHorse and choose candidate reference genomes
##    (all genomes in DarkHorse output)
bash ${scripts_dir}/run_darkhorse.sh --hit-table-fp ${hit_table_fp} \
                                     --darkhorse-config-fp ${darkhorse_config_fp} \
                                     --darkhose-install-dir ${darkhorse_install_dir} \
                                     --species_faa_fp ${species_faa_fp} \
                                     --working_dir ${working_dir} \
                                     --verbose ${verbose} \
                                     --lpi_upper ${lpi_upper} \
                                     --lpi_lower ${lpi_lower} \
                                     --parse_hgts false \
                                     --scripts_dir ${scripts_dir} \
                                     --output_fp ${output_fp}
selected_genomes=`
  cat ${working_dir}/darkhorse/calcs_*/*_smry | \
  sed -n '1!p' | \
  cut -f13 | \
  sort | \
  uniq`

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
                                       ${species_genbank_fp} \
                                       ${species_model_fp} \
                                       ${species_faa_fp} \
                                       ${ref_species_faa_fp} \
                                       ${gene_tree_dir} \
                                       ${gene_msa_dir} \
                                       ${hit_table_fp} \
                                       ${phylonet_install_dir} \
                                       ${jane_install_dir} \
                                       ${trex_install_dir} \
                                       ${verbose} \
                                       "${init_command}" \
                                       ${threads} \
                                       "${bash_config}" \
                                       ${qsub_env}
