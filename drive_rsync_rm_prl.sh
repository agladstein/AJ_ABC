#!/usr/bin/env bash

module load irods

set -e # quits at first error
set -x # print more debugging in stdout

OUT_PATH=$1

IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/8bfef8b8f3f7c43fd99832a323ef7130f98571bb/atmo_instance_ip.txt)

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

BUCKETS=$(find ${OUT_PATH} -maxdepth 1 -type d | tail -n +2 | rev | cut -d "/" -f1 | rev)
#BUCKETS=$(find ${OUT_PATH} -maxdepth 1 -type d | tail -n +2 | rev | cut -d "/" -f1 | rev | head -1)

for b in ${BUCKETS}; do
    # remove completed jobs
    CONTENT=$(find ${OUT_PATH}/$b | wc -l)
    SWITCH='keep'
    for q in ${JOBS}; do
        # if the job is not currently running
        if [ "$b" != "$q" ]; then
            # remove intermediate directories
            if [ -d "${OUT_PATH}/$b/sim_data_AJ_M1" ] || [ -d "${OUT_PATH}/$b/sim_data_AJ_M2" ] || [ -d "${OUT_PATH}/$b/sim_data_AJ_M3" ]; then
                echo "rm -r ${OUT_PATH}/$b/sim_data_AJ_M*"
                rm -r ${OUT_PATH}/$b/sim_data_AJ_M*
            fi
            if [ -d "${OUT_PATH}/$b/germline_out_AJ_M1" ] || [ -d "${OUT_PATH}/$b/germline_out_AJ_M2" ] || [ -d "${OUT_PATH}/$b/germline_out_AJ_M3" ]; then
                echo "rm -r ${OUT_PATH}/$b/germline_out_AJ_M*"
                rm -r ${OUT_PATH}/$b/germline_out_AJ_M*
            fi
            # if there aren't any contents of the output dirs
            if [ "${CONTENT}" -le 3 ]; then
                echo $b != $q and ${CONTENT} -le 3
                SWITCH='remove'
                continue
            fi
        fi
    done
    if [ ${SWITCH} == 'remove' ]; then
        echo "rm -r ${OUT_PATH}/$b"
        rm -r ${OUT_PATH}/$b
        continue
    fi

    MODEL=$(echo ${OUT_PATH}/${b}/results_sims_AJ_M* | rev | cut -d "_" -f1 | rev)

    # rsync to atmosphere
    ATMO_SIM_PATH=/vol_c/results_macsSwig_AJmodels_instant/sim_values_AJ_${MODEL}/${b}
    ATMO_RESULTS_PATH=/vol_c/results_macsSwig_AJmodels_instant/results_sims_AJ_${MODEL}/${b}

    ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_SIM_PATH}
    echo rsyncing ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/ ${ATMO_SIM_PATH}/
    rsync --remove-source-files -avz ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_SIM_PATH}/

    ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_RESULTS_PATH}
    echo rsyncing ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ ${ATMO_RESULTS_PATH}/
    rsync --remove-source-files -avz ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_RESULTS_PATH}/
done

exit


## Do the rest of the backing up from Atmosphere in different script.

# backup to google drive
DRIVE_SIM_PATH=backup_macsSwig_AJmodels_instant/sim_values_AJ_${MODEL}/${b}
DRIVE_RESULTS_PATH=backup_macsSwig_AJmodels_instant/results_sims_AJ_${MODEL}/${b}


echo 'google driving ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/ to '${DRIVE_SIM_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_SIM_PATH})" ]]; then
    echo ${DRIVE_SIM_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_SIM_PATH}
else
    echo ${DRIVE_SIM_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_SIM_PATH}/ ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/


echo 'google driving ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ to '${DRIVE_RESULTS_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_RESULTS_PATH})" ]]; then
    echo ${DRIVE_RESULTS_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_RESULTS_PATH}
else
    echo ${DRIVE_RESULTS_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_RESULTS_PATH}/ ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/



# backup to Data Store
IPLANT_SIM_PATH=/iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/sim_values_AJ_${MODEL}/${b}
IPLANT_RESULTS_PATH=/iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/results_sims_AJ_${MODEL}/${b}

echo 'iroding ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/'
imkdir -p ${IPLANT_SIM_PATH}
iput -K -f ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/ ${IPLANT_SIM_PATH}

echo 'iroding ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/'
imkdir -p ${IPLANT_RESULTS_PATH}
iput -K -f ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ ${IPLANT_RESULTS_PATH}