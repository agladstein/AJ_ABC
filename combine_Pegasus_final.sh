#!/usr/bin/env bash

VERSION=$1
CHR=$2

RESULTS_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}
ABC_DIR=/vol_c/ABC_AJmodels_${VERSION}

OSG_RESULTS=$(find ${RESULTS_DIR}/OSG -maxdepth 2 -type f | grep "chr${CHR}_")
CHTC_RESULTS=$(find ${RESULTS_DIR}/CHTC -maxdepth 2 -type f | grep "chr${CHR}_")

if [ -e ${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt ] ; then
    echo "${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt already exists"
    echo "rm ${ABC_DIR}/input_ABC_OSG_CHTC_${CHR}.txt"
    rm ${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt
fi
    
for RESULT in ${OSG_RESULTS} ${CHTC_RESULTS};
do
    echo ${RESULT}
    if [ -e ${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt ] ; then
	tail -n +2 ${RESULT} >>${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt
    else
	cat ${RESULT} >${ABC_DIR}/input_ABC_OSG_CHTC_chr${CHR}.txt
    fi
done
