#!/usr/bin/env bash

module load irods

set -e # quits at first error
set -x # print more debugging in stdout

#if [[ $HOSTNAME == service* ]]; then
#    module load irods
#fi

MODEL=$1
OUT_PATH=$2
SIM=$3
RESULTS=$4
PBS_ID=$5

IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/8bfef8b8f3f7c43fd99832a323ef7130f98571bb/atmo_instance_ip.txt)


# rsync to atmosphere
ATMO_SIM_PATH=/vol_c/results_macsSwig_AJmodels_instant/sim_values_AJ_M${MODEL}/${PBS_ID}
ATMO_RESULTS_PATH=/vol_c/results_macsSwig_AJmodels_instant/results_sims_AJ_M${MODEL}/${PBS_ID}

ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_SIM_PATH} # test
echo 'rsyncing ${OUT_PATH}/sim_values_AJ_M${MODEL}/'
rsync -a ${OUT_PATH}/sim_values_AJ_M${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_SIM_PATH}/

ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_RESULTS_PATH}
echo 'rsyncing ${OUT_PATH}/results_sims_AJ_M${MODEL}/'
rsync -a ${OUT_PATH}/results_sims_AJ_M${MODEL}/ agladstein@${IP_ADDRESS}:${ATMO_RESULTS_PATH}/


# backup to google drive
DRIVE_SIM_PATH=backup_macsSwig_AJmodels_instant/sim_values_AJ_M${MODEL}/$PBS_ID
DRIVE_RESULTS_PATH=backup_macsSwig_AJmodels_instant/results_sims_AJ_M${MODEL}/$PBS_ID

echo 'google driving ${OUT_PATH}/sim_values_AJ_M${MODEL}/ to '${DRIVE_SIM_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_SIM_PATH})" ]]; then
    echo ${DRIVE_SIM_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_SIM_PATH}
else
    echo ${DRIVE_SIM_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_SIM_PATH}/ ${OUT_PATH}/sim_values_AJ_M${MODEL}/

echo 'google driving ${OUT_PATH}/results_sims_AJ_M${MODEL}/ to '${DRIVE_RESULTS_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_RESULTS_PATH})" ]]; then
    echo ${DRIVE_RESULTS_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_RESULTS_PATH}
else
    echo ${DRIVE_RESULTS_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_RESULTS_PATH}/ ${OUT_PATH}/results_sims_AJ_M${MODEL}/


# backup to Data Store
echo 'iroding ${OUT_PATH}/sim_values_AJ_M${MODEL}/'
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/sim_values_AJ_M${MODEL}/${PBS_ID}
iput -K -f ${OUT_PATH}/sim_values_AJ_M${MODEL}/ /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/sim_values_AJ_M${MODEL}/${PBS_ID}

echo 'iroding ${OUT_PATH}/results_sims_AJ_M${MODEL}/'
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/results_sims_AJ_M${MODEL}/${PBS_ID}
iput -K -f ${OUT_PATH}/results_sims_AJ_M${MODEL}/ /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/results_sims_AJ_M${MODEL}/${PBS_ID}

echo 'rm file ' ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM}
rm ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM}

echo 'rm file ' ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS}
rm ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS}

