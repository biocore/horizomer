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
# usage: bash launch_software.sh working_dir scripts_dir species_tree_fp \
#             species_genome_fp species_model_fp query_species_coding_seqs_fp \
#             ref_species_coding_seqs_fp gene_tree_dir gene_msa_dir \
#             phylonet_install_dir jane_install_dir trex_install_dir \
#             verbose_str

# working dir
working_dir=$(readlink -m $1)
# scripts dir
scripts_dir=$(readlink -m $2)
# species tree in Newick format
species_tree_fp=$3
# species raw genome in FASTA format
species_genome_fp=$4
# species HMM model (produced by GeneMarkS)
species_model_fp=$5
# query species protein coding sequences in FASTA format 
query_species_coding_seqs_fp=$6
# reference species protein coding sequences in FASTA format (DB folder from ALF)
ref_species_coding_seqs_fp=$7
# gene trees in Newick format (GeneTrees folder from ALF)
gene_tree_dir=$8
# gene multiple sequence alignment dir
gene_msa_dir=$9
# DIAMOND alignments for query genome
diamond_tabular_query=${10}
# PhyloNet install dir
phylonet_install_dir=${11}
# Jane 4 install dir
jane_install_dir=${12}
# T-REX install dir
trex_install_dir=${13}
# Verbose string 'true' or 'false'
verbose=${14}
# Initial command that precedes call to software
# (example choosing virtualenv to workon)
init_command="${15}"
# Number of threads
threads=${16}
# nr FASTA file
nr_fp=${17}
# directory that contains names.dmp and nodes.dmp
nr_tax_dp=${18}
# nr DIAMOND database
diamond_db_nr=${19}
# Bash config file path (if None, default ~/.bash_profile)
bash_config="${20}"
# Launch on qsub cluster environment (true or false, if None, defaults to true)
qsub_env=${21}
# DarkHorse config file
darkhorse_config_fp=${22}
# DarkHorse install directory
darkhorse_install_dp=${23}
# HGTector config file
hgtector_config_fp=${24}
# HGTector install directory
hgtector_install_dp=${25}
# GI-to-TaxID translation table
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip
gi_to_taxid_fp=${26}
# qsub params
qsub="-q route -m abe -M jenya.kopylov@gmail.com -l nodes=1:ppn=${threads} -l walltime=230:00:00 -l pmem=10gb -l mem=20gb"

if [ "$verbose" == "true" ]
then
    echo "working dir: $working_dir"
    echo "scripts_dir: $scripts_dir"
    echo "species tree: $species_tree_fp"
    echo "species genome: $species_genome_fp"
    echo "HMM model: $species_model_fp"
    echo "query genomes: $query_species_coding_seqs_fp"
    echo "species proteome: $ref_species_coding_seqs_fp"
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
    echo "darkhorse_config_fp: $darkhorse_config_fp"
    echo "darkhorse_install_dp: $darkhorse_install_dp"
    echo "hgtector_config_fp: $hgtector_config_fp"
    echo "hgtector_install_dp: $hgtector_install_dp"
fi

if [ "${bash_config}" == "None" ]
then
    bash_config="~/.bash_profile"
fi
if [ "${qsub_env}" == "None" ]
then
    qsub_env="true"
fi

base_input_file_nwk="input_tree_nwk"
base_input_file_nex="input_tree_nex"
base_output_file="results_prog"
hgt_summary=$working_dir/"hgt_summary"
input_file_nwk=$working_dir/$base_input_file_nwk
input_file_nex=$working_dir/$base_input_file_nex
output_file=$working_dir/$base_output_file
stderr=$working_dir/"stderr"
stdout=$working_dir/"stdout"

mkdir -p "${working_dir}"
if [ "${init_command}" == "None" ]
then
    init_command="sleep 1"
fi

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
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_trex; sleep 2
else
    echo "${cmd}"
fi

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
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_ranger; sleep 2
else
    echo "${cmd}"
fi                                    

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
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_riata; sleep 2
else
    echo "${cmd}"
fi 

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
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_jane4; sleep 2
else
    echo "${cmd}"
fi

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
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_consel; sleep 2
else
    echo "${cmd}"
fi

## run GeneMark
cmd="${init_command}; \
      bash ${scripts_dir}/run_genemark.sh ${species_model_fp} \
                                          ${output_file}.gm.txt \
                                          ${species_genome_fp} \
                                          ${stdout}.gm.txt \
                                          ${stderr}.gm.txt \
                                          ${working_dir}"
if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_genemark; sleep 2
else
    echo "${cmd}"
fi

## run DarkHorse
cmd="${init_command}; \
      bash ${scripts_dir}/run_darkhorse.sh ${nr_fp} \
                                           ${diamond_db_nr} \
                                           ${diamond_tabular_query} \
                                           ${darkhorse_config_fp} \
                                           ${darkhorse_install_dp} \
                                           ${query_species_coding_seqs_fp} \
                                           ${working_dir} \
                                           ${threads}"

if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_darkhorse; sleep 2
else
    echo "${cmd}"
fi

## run HGTector
cmd="${init_command}; \
      bash ${scripts_dir}/run_hgtector.sh ${nr_fp} \
                                          ${diamond_db_nr} \
                                          ${diamond_tabular_query} \
                                          ${hgtector_config_fp} \
                                          ${hgtector_install_dp} \
                                          ${query_species_coding_seqs_fp} \
                                          ${working_dir} \
                                          ${threads} \
                                          ${tax_id} \
                                          ${nr_tax_dp} \
                                          ${gi_to_taxid_fp}"

if [ "${qsub_env}" == "true" ]
then
    echo "source ${bash_config}; \
          ${cmd}" | qsub $qsub -N run_hgtector; sleep 2
else
    echo "${cmd}"
fi

