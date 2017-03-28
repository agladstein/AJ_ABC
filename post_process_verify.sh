#!/usr/bin/env bash

SIM_PATH=$1
VERIFY_PATH=$2
MODEL=$3

for OUT in results_sims_AJ_M sim_values_AJ_M
do

    n=$(ls ${SIM_PATH}/${OUT}${MODEL}/ | wc -l)
    HEADER_FILE=$(ls ${SIM_PATH}/${OUT}${MODEL} | head -1)
    echo $n
    echo ${SIM_PATH}/${OUT}${MODEL}/${HEADER_FILE}

    echo "sim" >${VERIFY_PATH}/tmp1
    head -1 ${SIM_PATH}/${OUT}${MODEL}/${HEADER_FILE} >${VERIFY_PATH}/tmp2
    paste ${VERIFY_PATH}/tmp1 ${VERIFY_PATH}/tmp2 >${VERIFY_PATH}/${OUT}${MODEL}_${n}.txt

    rm ${VERIFY_PATH}/tmp1; rm ${VERIFY_PATH}/tmp2

    for f in ${SIM_PATH}/${OUT}${MODEL}/*; do tail -1 $f >>${VERIFY_PATH}/stats; done

    for f in ${SIM_PATH}/${OUT}${MODEL}/*; do echo $f >>${VERIFY_PATH}/name; done

    paste ${VERIFY_PATH}/name ${VERIFY_PATH}/stats >>${VERIFY_PATH}/${OUT}${MODEL}_${n}.txt

    rm ${VERIFY_PATH}/stats; rm ${VERIFY_PATH}/name

done