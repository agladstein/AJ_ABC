#!/usr/bin/env bash

set -e # quits at first error
set -x # print more debugging in stdout

MODEL=$1
OUT_PATH=$2
SIM=$3
RESULTS=$4
PBS_ID=$5

IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/8bfef8b8f3f7c43fd99832a323ef7130f98571bb/atmo_instance_ip.txt)


# backup to google drive
echo 'google driving ' ${SIM}
if [[ $(~/bin/drive ls backup_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/$PBS_ID | grep -q "cannot be found remotely") -eq 0 ]];then
    ~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination backup_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/ ${OUT_PATH}/${SIM}
else
    ~/bin/drive new -folder backup_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}
    ~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination backup_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/ ${OUT_PATH}/${SIM}
fi

echo 'google driving ' ${RESULTS}
if [[ $(~/bin/drive ls backup_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/$PBS_ID | grep -q "cannot be found remotely") -eq 0 ]];then
    ~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination backup_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}/ ${OUT_PATH}/${RESULTS}
else
    ~/bin/drive new -folder backup_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}
    ~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination backup_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}/ ${OUT_PATH}/${RESULTS}
fi


# backup to Data Store
echo 'iroding ' ${SIM}
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}
iput -K ${OUT_PATH}/${SIM} /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}

echo 'iroding ' ${RESULTS}
imkdir -p /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}
iput -K ${OUT_PATH}/${RESULTS} /iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}


# rsync to atmosphere
ssh agladstein@${IP_ADDRESS} mkdir -p /vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID} # test
echo 'rsyncing ' ${SIM}
rsync -a ${OUT_PATH}/${SIM} agladstein@${IP_ADDRESS}:/vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/${SIM}

ssh agladstein@${IP_ADDRESS} mkdir -p /vol_c/results_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}
echo 'rsyncing ' ${RESULTS}
rsync -a ${OUT_PATH}/${RESULTS} agladstein@${IP_ADDRESS}:/vol_c/results_macsSwig_AJmodels_mfloat/results_sims_AJ_M${MODEL}/${PBS_ID}/${RESULTS}

# check that the files are the same on hpc and atmosphere
echo 'atmo md5sum ' ${SIM} $(ssh agladstein@${IP_ADDRESS} md5sum /vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/${SIM} | cut -d " " -f1)
echo 'hpc md5sum ' ${SIM} $(md5sum ${OUT_PATH}/${SIM} | cut -d " " -f1)

if [ $(ssh agladstein@${IP_ADDRESS} md5sum /vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/${SIM} | cut -d " " -f1) == $(md5sum ${OUT_PATH}/${SIM} | cut -d " " -f1) ] ; then
    echo 'rm file ' ${OUT_PATH}/${SIM}
    rm ${OUT_PATH}/${SIM}
fi

echo 'atmo md5sum ' ${RESULTS} $(ssh agladstein@${IP_ADDRESS} md5sum /vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/${RESULTS} | cut -d " " -f1)
echo 'hpc md5sum ' ${RESULTS} $(md5sum ${OUT_PATH}/${RESULTS} | cut -d " " -f1)

if [ $(ssh agladstein@${IP_ADDRESS} md5sum /vol_c/results_macsSwig_AJmodels_mfloat/sim_values_AJ_M${MODEL}/${PBS_ID}/${RESULTS} | cut -d " " -f1) == $(md5sum ${OUT_PATH}/${RESULTS} | cut -d " " -f1) ] ; then
    echo 'rm file ' ${OUT_PATH}/${RESULTS}
    rm ${OUT_PATH}/${RESULTS}
fi