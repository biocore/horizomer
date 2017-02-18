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

# note: OrthoFinder supports Python 2 only. Previously, a special version of
# OrthoFinder 0.4 which supports Python 3 was obtained from the developers.
# It is outdated however. In order to run up-to-date versions of OrthoFinder,
# a separate Python 2 conda environment is created, and called from this Bash
# script, instead of directly called from the Python scripts.

# note: OrthoFinder requires BLAST+, and it has an option that allows the user
# to provide pre-computed BLAST results. Therefore, there is room for hacking
# this process by providing DIAMOND results. It is to be further explored.

set -eu
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

$verbose && echo "Running OrthoFinder .."
[ -z $CONDA_PATH_BACKUP ] || CONDA_PATH_BACKUP=
source activate ${py2_conda_env}

# command
cmd="orthofinder -f ${input_faa_dir} -t ${threads}"
$verbose && echo "Command:"$'\n'"  $cmd"

# run OrthoFinder and record time
TIMEFORMAT='%U %R'
TIME="$( time ($cmd 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $TIME | awk '{print $1}')
wall_time=$(echo $TIME | awk '{print $2}')
echo "Total user time OrthoFinder: ${user_time}" >> $stderr
echo "Total wall time OrthoFinder: ${wall_time}" >> $stderr

# transfer output files
outdir=${working_dir}/orthofinder
[ -d $outdir ] && rm -rf $outdir
mv ${input_faa_dir}/Results_* $outdir

source deactivate
$verbose && echo "Done"
