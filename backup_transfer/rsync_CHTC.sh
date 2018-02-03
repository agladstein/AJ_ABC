#!/usr/bin/env bash

# find each CHTC version directory with outputs (workflows that have completed).
# extract the workflow id
# create a directory on Atmo called the workflow id
# rsync the final_results from OSG to Atmo

VERSION=$1
CHTC_DIR=/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_$VERSION
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/CHTC

HOST_NAME=chtc

echo "#########################"
echo "ATMO ${ATMO_DIR}"
echo "CHTC ${CHTC_DIR}"

ssh ${HOST_NAME} find ${CHTC_DIR}*/outputs -maxdepth 1 -type f -name "final_results.txt"| cut -d "/" -f7 | cut -d "_" -f4- | xargs -n 1 -I % echo "mkdir -p ${ATMO_DIR}/macsswig_simsaj_$VERSION_% && rsync -avz agladstein@${HOST_NAME}:${OSG_DIR}_%/outputs/final_results.txt ${ATMO_DIR}/macsswig_simsaj_$VERSION_%/" | bash




