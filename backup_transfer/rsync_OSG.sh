#!/usr/bin/env bash

# find each OSG version directory with outputs (workflows that have completed).
# extract the workflow id
# create a directory on Atmo called the workflow id
# rsync the final_results from OSG to Atmo

VERSION=$1
OSG_DIR=/local-scratch/agladstein/workflows/macsswig_simsaj_$VERSION
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/OSG

HOST_NAME=login02.osgconnect.net

echo "#########################"
echo "ATMO ${ATMO_DIR}"
echo "OSG ${OSG_DIR}"

ssh agladstein@${HOST_NAME} find ${OSG_DIR}*/ -maxdepth 1 -type d -name "outputs" | rev | cut -d "_" -f1 | rev | cut -d "/" -f1 | xargs -n 1 -I % echo "mkdir -p ${ATMO_DIR}/% && rsync -avz agladstein@${HOST_NAME}:${OSG_DIR}_%/outputs/final_results.txt ${ATMO_DIR}/%/" | bash




