#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# purpose: launch HGT software on testing data sets, reformat input trees to
#          follow input format for each tool and parse output files to report
#          standardized statistics (number of HGTs, donors, recipients, gains
#          and losses)
#

# working dir
working_dir=$(readlink -m $1)
# scripts dir
scripts_dir=$(readlink -m $2)
# species tree in Newick format
species_tree_fp=$3
# reference genome in FASTA format (for compositional methods)
species_genome_fp=$4
# species HMM model (produced by GeneMarkS)
species_model_fp=$5
# query species protein coding sequences in FASTA format 
query_species_coding_seqs_fp=$6
# gene trees in Newick format (GeneTrees folder from ALF)
gene_tree_dir=$7
# gene multiple sequence alignment dir
gene_msa_dir=$8
# DIAMOND alignments for query genome
diamond_tabular_query_fp=$9
# PhyloNet install dir
phylonet_install_dir=${10}
# Jane 4 install dir
jane_install_dir=${11}
# T-REX install dir
trex_install_dir=${12}
# Verbose string 'true' or 'false'
verbose=${13}
# Initial command that precedes call to software
# (example choosing virtualenv to workon)
init_command="${14}"
# Number of threads
threads=${15}
# nr FASTA file
nr_fp=${16}
# directory that contains names.dmp and nodes.dmp
nr_tax_dp=${17}
# nr DIAMOND database
diamond_db_nr=${18}
# Bash config file path (if None, default ~/.bash_profile)
bash_config="${19}"
# Launch on qsub cluster environment (true or false, if None, defaults to true)
qsub_env=${20}
# HGTector config file
hgtector_config_fp=${21}
# HGTector install directory
hgtector_install_dp=${22}
# GI-to-TaxID translation table
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip
gi_to_taxid_fp=${23}
# EGID install directory
egid_install_dp=${24}
# HGT summary file basename
hgt_summary=${25}
# Qsub command
qsub="${26}"


if [ "$verbose" == "true" ]
then
    echo "working dir: $working_dir"
    echo "scripts_dir: $scripts_dir"
    echo "species tree: $species_tree_fp"
    echo "species genome: $species_genome_fp"
    echo "HMM model: $species_model_fp"
    echo "query genomes: $query_species_coding_seqs_fp"
    echo "gene trees: $gene_tree_dir"
    echo "gene MSAs: $gene_msa_dir"
    echo "DIAMOND alignments of query genome: ${diamond_query_tabular}"
    echo "PhyloNet install dir: $phylonet_install_dir"
    echo "Jane 4 install dir: $jane_install_dir"
    echo "T-REX install dir: $trex_install_dir"
    echo "verbose: $verbose"
    echo "initial command: $init_command"
    echo "DIAMOND nr: $diamond_db_nr"
    echo "nr : $nr_fp"
    echo "nr_tax_dp: $nr_tax_dp"
    echo "threads: $threads"
    echo "bash config file: $bash_config"
    echo "qsub_env: $qsub_env" 
    echo "hgtector_config_fp: $hgtector_config_fp"
    echo "hgtector_install_dp: $hgtector_install_dp"
    echo "egid_install_dp: $egid_install_dp"
fi

base_input_file_nwk="input_tree_nwk"
base_input_file_nex="input_tree_nex"
base_output_file="results_prog"
input_file_nwk=$working_dir/$base_input_file_nwk
input_file_nex=$working_dir/$base_input_file_nex
output_file=$working_dir/$base_output_file
stderr=$working_dir/"stderr"
stdout=$working_dir/"stdout"

## run T-REX
cmd="${init_command}; \
      bash ${scripts_dir}/run_trex.sh ${gene_tree_dir} \
                                      ${hgt_summary}.trex.txt \
                                      ${verbose}\
                                      ${stdout}.trex.txt \
                                      ${stderr}.trex.txt \
                                      ${scripts_dir} \
                                      ${species_tree_fp} \
                                      ${input_file_nwk}.trex.txt \
                                      ${trex_install_dir} \
                                      ${base_input_file_nwk}.trex.txt"
submit_job "${cmd}" trex

## run RANGER-DTL
cmd="${init_command}; \
      bash ${scripts_dir}/run_ranger.sh ${gene_tree_dir} \
                                        ${hgt_summary}.ranger.txt \
                                        ${verbose} \
                                        ${stdout}.ranger.txt \
                                        ${stderr}.ranger.txt \
                                        ${scripts_dir} \
                                        ${species_tree_fp} \
                                        ${input_file_nwk}.ranger.txt \
                                        ${output_file}.ranger.txt"
submit_job "${cmd}" ranger                                   

## run RIATA-HGT
cmd="${init_command}; \
      bash ${scripts_dir}/run_riatahgt.sh ${gene_tree_dir} \
                                          ${hgt_summary}.riatahgt.txt \
                                          $verbose \
                                          ${stdout}.riatahgt.txt \
                                          ${stderr}.riatahgt.txt \
                                          ${scripts_dir} \
                                          ${species_tree_fp} \
                                          ${input_file_nex}.riata.txt \
                                          ${output_file}.riatahgt.txt \
                                          ${phylonet_install_dir}"
submit_job "${cmd}" riatahgt

## run JANE 4
cmd="${init_command}; \
      bash ${scripts_dir}/run_jane4.sh ${gene_tree_dir} \
                                       ${hgt_summary}.jane4.txt \
                                       $verbose \
                                       ${stdout}.jane4.txt \
                                       ${stderr}.jane4.txt \
                                       ${scripts_dir} \
                                       ${species_tree_fp} \
                                       ${input_file_nex}.jane.txt \
                                       ${output_file}.jane4.txt \
                                       ${jane_install_dir}"
submit_job "${cmd}" jane4

## run CONSEL
cmd="${init_command}; \
      bash ${scripts_dir}/run_consel.sh ${gene_tree_dir} \
                                        ${hgt_summary}.consel.txt \
                                        $verbose \
                                        ${stdout}.consel.txt \
                                        ${stderr}.consel.txt \
                                        ${scripts_dir} \
                                        ${species_tree_fp} \
                                        ${input_file_nwk}.consel.txt \
                                        ${output_file}.consel.txt \
                                        ${gene_msa_dir} \
                                        ${working_dir} "
submit_job "${cmd}" consel

## run GeneMark
cmd="${init_command}; \
      bash ${scripts_dir}/run_genemark.sh ${species_model_fp} \
                                          ${output_file}.gm.txt \
                                          ${species_genome_fp} \
                                          ${stdout}.gm.txt \
                                          ${stderr}.gm.txt \
                                          ${working_dir}"
submit_job "${cmd}" genemark

## run HGTector
cmd="${init_command}; \
      bash ${scripts_dir}/run_hgtector.sh ${nr_fp} \
                                          ${diamond_db_nr} \
                                          ${diamond_tabular_query_fp} \
                                          ${hgtector_config_fp} \
                                          ${hgtector_install_dp} \
                                          ${query_species_coding_seqs_fp} \
                                          ${working_dir} \
                                          ${threads} \
                                          ${tax_id} \
                                          ${nr_tax_dp} \
                                          ${gi_to_taxid_fp} \
                                          ${scripts_dir} \
                                          ${hgt_summary}.hgtector.txt"
submit_job "${cmd}" hgtector

## run EGID
cmd="${init_command}; \
      bash ${scripts_dir}/run_egid.sh ${species_genome_fp} \
                                      ${scripts_dir} \
                                      ${working_dir} \
                                      ${egid_install_dir} \
                                      ${hgt_summary}.egid.txt"
submit_job "${cmd}" egid


