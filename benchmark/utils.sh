# launch job
function submit_job {
    cmd=${1}
    tool=${2}
    if [ "${qsub_env}" == "true" ]
    then
        echo "source ${bash_config}; \
             ${cmd}" | qsub $qsub -N run_$tool; sleep 2
    else
        echo "${cmd}"
    fi
}