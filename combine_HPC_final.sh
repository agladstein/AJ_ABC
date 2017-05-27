#!/usr/bin/env bash

OUT_PATH=$1

for bucket in $(ls ${OUT_PATH})
do
    PRE=$(echo $bucket | cut -d "_" -f1-3)
    BUCKET_FILE=${OUT_PATH}/$bucket
    FILE=${OUT_PATH}/../${PRE}_HPC.txt
    if [ -e ${FILE} ] ; then
        grep -v "Asc_NAF" ${BUCKET_FILE} >>${FILE}
    else
        cp ${BUCKET_FILE} ${FILE}
    fi
done
