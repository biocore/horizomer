#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run OrthoFinder to infer orthologous groups from protein-coding genes
# of selected genomes

# note: OrthoFinder supports Python 2 only. Therefore, a separate Python 2
# conda environment, as defined by variable $py2_conda_env, is activated in
# order to run the program.

# note: OrthoFinder requires BLAST+ for homology search. Meanwhile, it has an
# option that allows the user to provide pre-computed BLAST results. Therefore,
# I hacked this process by mocking BLAST results with DIAMOND results.

set -e  # no -u otherwise causes error when switching conda
source $(dirname "$0")/utils.sh
args=(
    working_dir
    input_faa_dir
    py2_conda_env
    threads
    stdout
    stderr
    verbose
)
get_args "$@"

$verbose && echo "Preparing files for OrthoFinder .."
ofdir=${working_dir}/orthofinder
mkdir -p $ofdir

# convert species names and sequence names into incremental numbers, as
# required by OrthoFinder
isp=-1  # index of current species
iseq=0  # index of current sequence
for sp in ${input_faa_dir}/*.faa
do
    ((isp++))
    echo $isp': '$(basename $sp) >> $ofdir/SpeciesIDs.txt
    iseq=0
    while read line
    do
      	if [[ $line == '>'* ]]
        then
            echo $isp'_'$iseq': '${line#>} >> $ofdir/SequenceIDs.txt
            echo '>'$isp'_'$iseq >> $ofdir/Species$isp.fa
            ((iseq++))
        else
            echo $line >> $ofdir/Species$isp.fa
        fi
    done < $sp
    diamond makedb --in $ofdir/Species$isp.fa \
                   --db $ofdir/Species$isp \
                   --threads $threads
                   --quiet
done
$verbose && echo "Done"

# run all-to-all DIAMOND searches to infer homology
$verbose && echo "Running DIAMOND .."
for i in $(seq 0 $isp)
do
    for j in $(seq 0 $isp)
    do
      	diamond blastp --query $ofdir/Species$i.fa \
                       --db $ofdir/Species$j \
                       --evalue 0.001 \
                       --out $ofdir/Blast$i'_'$j.txt \
                       --threads $threads
                       --quiet
    done
done
$verbose && echo "Done"

$verbose && echo "Running OrthoFinder .."
source activate ${py2_conda_env}

# command
cmd="orthofinder --blast $ofdir --threads $threads"
$verbose && echo "Command:"$'\n'"  $cmd"

# run OrthoFinder and record time
TIMEFORMAT='%U %R'
TIME="$( time ($cmd 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $TIME | awk '{print $1}')
wall_time=$(echo $TIME | awk '{print $2}')
echo "Total user time OrthoFinder: ${user_time}" >> $stderr
echo "Total wall time OrthoFinder: ${wall_time}" >> $stderr

source deactivate
$verbose && echo "Done"
