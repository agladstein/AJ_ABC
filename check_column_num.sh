#!/usr/bin/env bash

#Does not work right now.

#make file with number of desired columns

FILE=$1
COL=$2

NUM=$(wc -l ${FILE} | cut -d " " -f1)
for i in {1..${NUM}}
do
    if [ $(head -${i} ${FILE} | tail -1 | tr '\t' '\n' | wc -l) -eq '${COL}' ];then
        echo "hello"
#        head -${i} ${FILE} | tail -1 >>tmp1
    fi
done
#chmv tmp1 ${FILE}