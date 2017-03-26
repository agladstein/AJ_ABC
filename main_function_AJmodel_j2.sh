#!/usr/bin/env bash
SYSTEM=$1

if [ "$SYSTEM" == "ocelote" ]; then
    VIRTUAL_ENV=/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/macss_env
    PATH=/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/macss_env/bin:$PATH
    for i in {1..3}
    do
	QUE=standard JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=macss_env/bin/python MODEL=${i} j2 template.pbs.j2 >model${i}_${SYSTEM}_standard.pbs
	QUE=qualified JNUM=1-2500 NODE=1 CORE=1 MEM=6 PYTHONENV=macss_env/bin/python MODEL=${i} j2 template.pbs.j2 >model${i}_${SYSTEM}_qualified.pbs
	QUE=windfall JNUM=1-5000 NODE=1 CORE=1 MEM=6 PYTHONENV=macss_env/bin/python MODEL=${i} j2 template.pbs.j2 >model${i}_${SYSTEM}_windfall.pbs
    done
else
    VIRTUAL_ENV=/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/macss_env_ice
    PATH=/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/macss_env_ice/bin:$PATH
    for i in {1..3}
    do
	QUE=standard JNUM=1-1500 NODE=1 CORE=1 MEM=4 PYTHONENV=macss_env_ice/bin/python MODEL=${i} JOBTYPE=${SYSTEM}_only j2 template.pbs.j2 >model${i}_${SYSTEM}_standard.pbs
	QUE=qualified JNUM=1-1500 NODE=1 CORE=1 MEM=4 PYTHONENV=macss_env_ice/bin/python MODEL=${i} JOBTYPE=${SYSTEM}_only j2 template.pbs.j2 >model${i}_${SYSTEM}_qualified.pbs
	QUE=windfall JNUM=1-5000 NODE=1 CORE=1 MEM=4 PYTHONENV=macss_env_ice/bin/python MODEL=${i} JOBTYPE=${SYSTEM}_only j2 template.pbs.j2 >model${i}_${SYSTEM}_windfall.pbs
    done
fi







