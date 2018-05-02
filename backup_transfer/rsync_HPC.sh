#!/usr/bin/env bash

set -e # quits at first error

VERSION=$1
HPC_DIR=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_${VERSION}
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/HPC

HOST_NAME=sftp.hpc.arizona.edu

echo "#########################"
echo "ATMO ${ATMO_DIR}"
echo "HPC ${HPC_DIR}"

echo "rsync -azv --remove-source-files --exclude '*Posterior*' --exclude '*germline*' --exclude '*sim_data*' --exclude '*ped' --exclude '*map' --exclude '*match' --exclude '*log' --exclude 'macsargs*' agladstein@${HOST_NAME}:${HPC_DIR}/ ${ATMO_DIR}"
rsync -azv --remove-source-files --exclude '*Posterior*' --exclude '*germline*' --exclude '*sim_data*' --exclude '*ped' --exclude '*map' --exclude '*match' --exclude '*log' --exclude 'macsargs*' agladstein@${HOST_NAME}:${HPC_DIR}/ ${ATMO_DIR}


echo "Finished transfer"
