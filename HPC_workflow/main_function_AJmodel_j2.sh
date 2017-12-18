#!/usr/bin/env bash

set -e # quits at first error

SYSTEM=$1
OUT=$2
MODEL=$3

if [ "$SYSTEM" == "ocelote" ]; then
#    PYTHONMOD=python/2/2.7.11
#    module load ${PYTHONMOD}
#    PBS_ENV=/home/u15/agladstein/env/macss_env_ocelote_2.7.11
#    module load python/3/3.5.2
#    VIRTUAL_ENV=/12kx/gsfs2/home/u15/agladstein/ABC/macsSwig_AJmodels/j2_env_ocelote
#    PATH=${VIRTUAL_ENV}/bin:$PATH
#	PYTHONMOD=${PYTHONMOD} QUE=standard JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=${PBS_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_standard.pbs
#	PYTHONMOD=${PYTHONMOD} QUE=qualified JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=${PBS_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_qualified.pbs
#	PYTHONMOD=${PYTHONMOD} QUE=windfall JNUM=1-2000 NODE=1 CORE=1 MEM=6 PYTHONENV=${PBS_ENV}/bin/python MODEL=${MODEL} OUT_PATH=${OUT} j2 template.pbs.j2 >model${MODEL}_${SYSTEM}_windfall.pbs
    echo "skipping jinja, pbs already made"
else
    PYTHONMOD=python/2.7.9
    module load ${PYTHONMOD}
    export MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles:/uaopt/modulefiles
    eval `/usr/bin/modulecmd bash load ${PYTHONMOD}`
    VIRTUAL_ENV=/home/u15/agladstein/env/macss_env_ICE_2.7.9
    PATH=${VIRTUAL_ENV}/bin:$PATH
    if [ "$SYSTEM" == "smp" ]; then
        ICE_MEM=8
    else
        ICE_MEM=4
    fi
	PYTHONMOD=${PYTHONMOD} QUE=windfall JNUM=1-1000 TIME=1 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 macsargs_pbs.j2 >PBS/macsargs_model${MODEL}_${SYSTEM}_windfall.pbs
	PYTHONMOD=${PYTHONMOD} QUE=windfall JNUM=1-22 TIME=60 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 run-sim_pbs.j2 >PBS/run-sim_model${MODEL}_${SYSTEM}_windfall.pbs
	PYTHONMOD=${PYTHONMOD} QUE=windfall TIME=5 NODE=1 CORE=1 MEM=${ICE_MEM} PYTHONENV=${VIRTUAL_ENV}/bin/python MODEL=${MODEL} JOBTYPE=${SYSTEM}_only OUT_PATH=${OUT} j2 genome_stats_pbs.j2 >PBS/genome_stats_model${MODEL}_${SYSTEM}_windfall.pbs
fi







