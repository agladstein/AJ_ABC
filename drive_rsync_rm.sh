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

if [ ! -f ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} ] && [ ! -f ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} ]; then
    echo ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} 'and' ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} 'not found'
    exit 1
fi

# rsync to atmosphere
ATMO_SIM_PATH=/vol_c/results_macsSwig_AJmodels_rscale4Trel100/sim_values_AJ_M${MODEL}/${PBS_ID}
ATMO_RESULTS_PATH=/vol_c/results_macsSwig_AJmodels_rscale4Trel100/results_sims_AJ_M${MODEL}/${PBS_ID}

ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_SIM_PATH} # test
echo 'rsyncing ' ${SIM}
rsync -a ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} agladstein@${IP_ADDRESS}:${ATMO_SIM_PATH}/${SIM}

ssh agladstein@${IP_ADDRESS} mkdir -p ${ATMO_RESULTS_PATH}
echo 'rsyncing ' ${RESULTS}
rsync -a ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} agladstein@${IP_ADDRESS}:${ATMO_RESULTS_PATH}/${RESULTS}


# backup to google drive
DRIVE_SIM_PATH=backup_macsSwig_AJmodels_rscale4Trel100/sim_values_AJ_M${MODEL}/$PBS_ID
DRIVE_RESULTS_PATH=backup_macsSwig_AJmodels_rscale4Trel100/results_sims_AJ_M${MODEL}/$PBS_ID

echo 'google driving ' ${SIM} 'to '${DRIVE_SIM_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_SIM_PATH})" ]]; then
    echo ${DRIVE_SIM_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_SIM_PATH}
else
    echo ${DRIVE_SIM_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_SIM_PATH}/ ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM}

echo 'google driving ' ${RESULTS} 'to '${DRIVE_RESULTS_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_RESULTS_PATH})" ]]; then
    echo ${DRIVE_RESULTS_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_RESULTS_PATH}
else
    echo ${DRIVE_RESULTS_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_RESULTS_PATH}/ ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS}


# backup to Data Store
echo 'iroding ' ${SIM}
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_rscale4Trel100/sim_values_AJ_M${MODEL}/${PBS_ID}
iput -K -f ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_rscale4Trel100/sim_values_AJ_M${MODEL}/${PBS_ID}

echo 'iroding ' ${RESULTS}
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_rscale4Trel100/results_sims_AJ_M${MODEL}/${PBS_ID}
iput -K -f ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_rscale4Trel100/results_sims_AJ_M${MODEL}/${PBS_ID}



# check that the files are the same on hpc and atmosphere
echo 'atmo md5sum ' ${SIM} $(ssh agladstein@${IP_ADDRESS} md5sum ${ATMO_SIM_PATH}/${SIM} | cut -d " " -f1)
echo 'hpc md5sum ' ${SIM} $(md5sum ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} | cut -d " " -f1)

if [[ $(ssh agladstein@${IP_ADDRESS} md5sum ${ATMO_SIM_PATH}/${SIM} | cut -d " " -f1) == $(md5sum ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM} | cut -d " " -f1) ]] ; then
    echo 'rm file ' ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM}
    rm ${OUT_PATH}/sim_values_AJ_M${MODEL}/${SIM}
fi

echo 'atmo md5sum ' ${RESULTS} $(ssh agladstein@${IP_ADDRESS} md5sum ${ATMO_RESULTS_PATH}/${RESULTS} | cut -d " " -f1)
echo 'hpc md5sum ' ${RESULTS} $(md5sum ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} | cut -d " " -f1)

if [[ $(ssh agladstein@${IP_ADDRESS} md5sum ${ATMO_RESULTS_PATH}/${RESULTS} | cut -d " " -f1) == $(md5sum ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS} | cut -d " " -f1) ]] ; then
    echo 'rm file ' ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS}
    rm ${OUT_PATH}/results_sims_AJ_M${MODEL}/${RESULTS}
fi