#!/usr/bin/env bash

set -e # quits at first error

VERSION=$1
HPC_DIR=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_${VERSION}
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/HPC


IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/669a0e602306776ffa8d8be33e63574dfa2d1766/atmo_instance_ip.txt)

# check for currently running jobs
IP_service1=128.196.131.51
IP_login2=150.135.165.106
if [[ $HOSTNAME == service* ]]; then
    ICE_JOBS=$(/usr/local/bin/qstat_local | grep "agladstein" | cut -d "[" -f1)
    OCELOTE_JOBS=$(ssh agladstein@${IP_login2} /cm/shared/apps/pbspro/current/bin/qstat | grep "agladstein" | cut -d "[" -f1)
else
    ICE_JOBS=$(ssh agladstein@${IP_service1} /usr/local/bin/qstat_local | grep "agladstein" | cut -d "[" -f1)
    OCELOTE_JOBS=$(/cm/shared/apps/pbspro/current/bin/qstat | grep "agladstein" | cut -d "[" -f1)
fi
JOBS=$OCELOTE_JOBS' '$ICE_JOBS

BUCKETS=$(find ${HPC_DIR} -maxdepth 1 -type d | tail -n +2 | rev | cut -d "/" -f1 | rev)

for BUCKET in ${BUCKETS}; do
    # remove completed jobs
    CONTENT=$(find ${HPC_DIR}/${BUCKER} | wc -l)
    SWITCH='keep'
    for q in ${JOBS}; do
        # if the job is not currently running
        if [ "${BUCKER}" != "$q" ]; then
            # remove intermediate directories
            if [ -d "${HPC_DIR}/${BUCKER}/sim_data_AJ_M1" ] || [ -d "${HPC_DIR}/${BUCKER}/sim_data_AJ_M2" ] || [ -d "${HPC_DIR}/${BUCKER}/sim_data_AJ_M3" ]; then
                echo "rm -r ${HPC_DIR}/${BUCKER}/sim_data_AJ_M*"
                rm -r ${HPC_DIR}/${BUCKER}/sim_data_AJ_M*
            fi
            if [ -d "${HPC_DIR}/${BUCKER}/germline_out_AJ_M1" ] || [ -d "${HPC_DIR}/${BUCKER}/germline_out_AJ_M2" ] || [ -d "${HPC_DIR}/${BUCKER}/germline_out_AJ_M3" ]; then
                echo "rm -r ${HPC_DIR}/${BUCKER}/germline_out_AJ_M*"
                rm -r ${HPC_DIR}/${BUCKER}/germline_out_AJ_M*
            fi
            # if there aren't any contents of the output dirs
            if [ "${CONTENT}" -le 3 ]; then
                echo ${BUCKER} != $q and ${CONTENT} -le 3
                SWITCH='remove'
                continue
            fi
        fi
    done
    if [ ${SWITCH} == 'remove' ]; then
        echo "rm -r ${HPC_DIR}/${BUCKER}"
        rm -r ${HPC_DIR}/${BUCKER}
        continue
    fi

    MODEL=$(echo ${HPC_DIR}/${BUCKET}/results_sims_AJ_M* | rev | cut -d "_" -f1 | rev)

    # rsync to atmosphere
    ATMO_SIM_PATH=/vol_c/results_macsSwig_AJmodels_instant/sim_values_AJ_${MODEL}/${BUCKET}
    ATMO_RESULTS_PATH=/vol_c/results_macsSwig_AJmodels_instant/results_sims_AJ_${MODEL}/${BUCKET}

    ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_SIM_PATH}
    echo rsyncing ${HPC_DIR}/${BUCKET}/sim_values_AJ_${MODEL}/ ${ATMO_SIM_PATH}/
    rsync --remove-source-files -avz ${HPC_DIR}/${BUCKET}/sim_values_AJ_${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_SIM_PATH}/

    ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_RESULTS_PATH}
    echo rsyncing ${HPC_DIR}/${BUCKET}/results_sims_AJ_${MODEL}/ ${ATMO_RESULTS_PATH}/
    rsync --remove-source-files -avz ${HPC_DIR}/${BUCKET}/results_sims_AJ_${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_RESULTS_PATH}/
done


