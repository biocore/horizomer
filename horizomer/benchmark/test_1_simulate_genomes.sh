#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# purpose: simulate genomes with ALF
# usage  : bash test_1_simulate_genomes.sh root_genome.fasta custom_tree.nwk scripts_dir working_dir

root_genome_fp=$(readlink -m $1)
custom_tree_fp=$(readlink -m $2)
scripts_dir=$(readlink -m $3)
working_dir=$(readlink -m $4)
alf_params="alf_params.txt"
stderr="stderr.log"
stdout="stdout.log"

lgt_rate=0.005
orth_rep_a=(1 0.5)
gc_cont_am_a=("False" "True")
gene_loss_rate_a=(0 0.005)
gene_dup_rate_a=(0 0.0006)
if [ ! -d "${working_dir}" ]; then
    mkdir $working_dir
fi

i=0
echo "Begin simulation .."
for orth_rep in "${orth_rep_a[@]}"
do
    echo -e "\tlgt rate: ${lgt_rate}"
    echo -e "\torth_rep: ${orth_rep}"
    for gc_cont_am in "${gc_cont_am_a[@]}"
    do
        echo -e "\tgc_content: ${gc_cont_am}"
        for gene_loss_rate in "${gene_loss_rate_a[@]}"
        do
            echo -e "\tgene loss rate: ${gene_loss_rate}"
            for gene_dup_rate in "${gene_dup_rate_a[@]}"
            do
                echo -e "\tgene duplication rate: ${gene_dup_rate}"
                echo -e "\toutput directory: ${working_dir}/params_$i"
                if [ ! -d "${working_dir}/params_$i" ]; then
                    mkdir ${working_dir}/"params_$i"
                fi
                python $scripts_dir/create_alf_params.py ${root_genome_fp} \
                                                         ${custom_tree_fp} \
                                                         ${working_dir}/"params_$i" \
                                                         ${alf_params} \
                                                         ${lgt_rate} \
                                                         ${orth_rep} \
                                                         ${gc_cont_am} \
                                                         ${gene_loss_rate} \
                                                         ${gene_dup_rate} \
                                                         "params_$i"

                printf "p(gene loss)\tp(gene duplication)\tp(orthologous gene replacement)\tGC content amelioration\n" > ${working_dir}/"params_$i"/"parameters_summary.txt"
                printf "${gene_loss_rate}\t${gene_dup_rate}\t${orth_rep}\t${gc_cont_am}\n" >> ${working_dir}/"params_$i"/"parameters_summary.txt"
                # launch ALF
                echo -e "\tRunning ALF .."
                (cd ${working_dir}/"params_$i"; alfsim "./alf_params.txt" 1>$stdout 2>$stderr)

                # format the ALF genes tree (Newick) to replace '/' with '_' and
                # remove the "[&&NHX:D=N]" tags
                echo -e "Cleaning Newick files .."
                for file in ${working_dir}/"params_$i"/"params_$i"/GeneTrees/*.nwk;
                do
                    sed -i "s/\//\_/g" $file
                    sed -i "s/\[&&NHX:D=N\]//g" $file
                    # remove empty lines
                    sed -i "/^$/d" $file
                done
                i=$((i+1))
            done
        done
    done
done
