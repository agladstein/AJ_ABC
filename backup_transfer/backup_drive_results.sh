#!/usr/bin/env bash

set -e # quits at first error

echo 'google driving ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/ to '${DRIVE_RESULTS_PATH}
if [[ -z "$(~/bin/drive ls ${DRIVE_RESULTS_PATH})" ]]; then
    echo ${DRIVE_RESULTS_PATH} 'is not made yet on google drive'
    ~/bin/drive new -folder ${DRIVE_RESULTS_PATH}
else
    echo ${DRIVE_RESULTS_PATH} 'is already made on google drive'
fi
~/bin/drive push -verbose -exclude-ops "delete,update" -no-prompt -destination ${DRIVE_RESULTS_PATH}/ ${OUT_PATH}/${b}/results_sims_AJ_${MODEL}/
