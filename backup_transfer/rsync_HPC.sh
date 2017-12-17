#!/usr/bin/env bash

set -e # quits at first error

VERSION=$1
HPC_DIR=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_genome/${VERSION}
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/HPC

HOST_NAME=sftp.hpc.arizona.edu
IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/669a0e602306776ffa8d8be33e63574dfa2d1766/atmo_instance_ip.txt)

echo "#########################"
echo "ATMO ${ATMO_DIR}"
echo "HPC ${HPC_DIR}"


rsync -avzn agladstein@${HOST_NAME}
