#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: run PhyloPhlAn to build phylogenetic tree of whole genomes

# prerequisites:
# one needs to install the development version of PhyloPhlAn:
#   hg clone https://bitbucket.org/nsegata/phylophlan
#   cd phylophlan
#   hg up dev
# one needs the follow programs installed and callable from $PATH:
#   usearch
#   muscle
#   FastTree
# one also needs biopython installed

set -eu
source $(dirname "$0")/utils.sh
args=(
    working_dir
    input_faa_dir
    scripts_dir
    phylophlan_install_dir
    threads
    stdout
    stderr
    verbose
)
get_args "$@"

if [ "$verbose" == "true" ]
then
    echo "Running PhyloPhlAn .."
fi

# set up
mkdir -p ${working_dir}/phylophlan
PWD=$(pwd)
cd ${working_dir}/phylophlan
mkdir -p input
mkdir -p output
mkdir -p temp
ln -s ${phylophlan_install_dir}/data
ln -s ${input_faa_dir} input/genomes

# command
cmd="${phylophlan_install_dir}/phylophlan.py -u input --nproc $threads --c_dat temp"
if [ "$verbose" == "true" ]
then
    echo "Command:"$'\n'"  $cmd"
fi

# run PhyloPhlAn and record time
TIMEFORMAT='%U %R'
TIME="$( time ($cmd 1>>$stdout 2>>$stderr) 2>&1)"
user_time=$(echo $TIME | awk '{print $1}')
wall_time=$(echo $TIME | awk '{print $2}')
echo "Total user time PhyloPhlAn: ${user_time}" >> $stderr
echo "Total wall time PhyloPhlAn: ${wall_time}" >> $stderr

# tear down
rm -rf temp
cd $PWD

# output genome tree will be:
# ${working_dir}/phylophlan/output/genomes.tree.nwk
