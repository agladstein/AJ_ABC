#!/usr/bin/env bash

VERSION=$1
MODEL=$2

RESULTS_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}
ABC_DIR=/vol_c/ABC_AJmodels_${VERSION}

find ${RESULTS_DIR}/HPC/model${MODEL} -type f | head -1 | xargs head -1 >${ABC_DIR}/input_ABC_HPC_${MODEL}.txt
find ${RESULTS_DIR}/HPC/model${MODEL} -type f | xargs -I % tail -n +2 % >>${ABC_DIR}/input_ABC_HPC_${MODEL}.txt
