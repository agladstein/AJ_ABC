#!/usr/bin/env bash

IN_PATH=$1
OUT_PATH=$2

for MODEL in {1..3}
do
    if [ "${MODEL}"=="1" ]; then
        HEADER=header_M${MODEL}_222.txt
    else
        HEADER=header_M${MODEL}.txt
    fi
    for BUCKET in $(ls ${IN_PATH}/sim_values_AJ_M${MODEL})
    do
        echo "python /vol_c/src/macsswig_simsaj/post_process.py ${IN_PATH} ${OUT_PATH} ${MODEL} ${BUCKET} ${HEADER} same &"
        python /vol_c/src/macsswig_simsaj/post_process.py ${IN_PATH} ${OUT_PATH} ${MODEL} ${BUCKET} ${HEADER} same &

        # if it starts to be too slow. do something else fancy to not make files that have already been made.
#        if [ ! -f ${OUT_PATH}/input_ABCtoolbox_M${MODEL}_${BUCKET}.txt ]; then
#            echo "python /vol_c/src/macsswig_simsaj/post_process.py ${IN_PATH} ${OUT_PATH} ${MODEL} ${BUCKET} ${HEADER} same &"
#	        python /vol_c/src/macsswig_simsaj/post_process.py ${IN_PATH} ${OUT_PATH} ${MODEL} ${BUCKET} ${HEADER} same &
#        else
#            echo "${OUT_PATH}/input_ABCtoolbox_M${MODEL}_${BUCKET}.txt is already created"
#        fi
    done
done
