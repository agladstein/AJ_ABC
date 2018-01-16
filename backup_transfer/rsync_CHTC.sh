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

ssh ${HOST_NAME} find ${CHTC_DIR}*/ -maxdepth 1 -type d -name "outputs" | rev | cut -d "_" -f1 | rev | cut -d "/" -f1 | xargs -n 1 -I % echo "mkdir -p ${ATMO_DIR}/% && rsync -navz ${HOST_NAME}:${CHTC_DIR}_%/outputs/final_results.txt ${ATMO_DIR}/%/" | bash




