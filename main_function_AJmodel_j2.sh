#!/usr/bin/env bash
SYSTEM=$1
OUT=$2
MODEL=$3

if [ "$SYSTEM" == "ocelote" ]; then
    VIRTUAL_ENV=/home/u15/agladstein/env/macss_env_ocelote_2.7.11
    PATH=${VIRTUAL_ENV}/bin:$PATH
	QUE=standard JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_standard.pbs
	QUE=qualified JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_qualified.pbs
	QUE=windfall JNUM=1-5000 NODE=1 CORE=1 MEM=6 PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_windfall.pbs

else
    VIRTUAL_ENV=/home/u15/agladstein/env/macss_env_ICE_2.7.9
    PATH=${VIRTUAL_ENV}/bin:$PATH
    if [ "$SYSTEM" == "smp" ]; then
        ICE_MEM=8
	    QUE=${SYSTEM}_qual JNUM=1-1500 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_qualified.pbs
    else
        ICE_MEM=4
    fi
    QUE=standard JNUM=1-1500 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_standard.pbs
	QUE=windfall JNUM=1-5000 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_windfall.pbs
fi







