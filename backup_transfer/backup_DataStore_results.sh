#!/usr/bin/env bash

set -e # quits at first error

OUT_PATH=$1

# backup to Data Store
IPLANT_SIM_PATH=/iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/sim_values_AJ_${MODEL}/${b}
IPLANT_RESULTS_PATH=/iplant/home/agladstein/AJmacs_data/macsSwig_AJmodels_instant/results_sims_AJ_${MODEL}/${b}

echo 'iroding ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/'
imkdir -p ${IPLANT_SIM_PATH}
iput -K -f ${OUT_PATH}/${b}/sim_values_AJ_${MODEL}/ ${IPLANT_SIM_PATH}

echo 'iroding ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/'
imkdir -p ${IPLANT_RESULTS_PATH}
iput -K -f ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ ${IPLANT_RESULTS_PATH}