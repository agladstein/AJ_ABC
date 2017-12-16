#!/usr/bin/env bash
cd /home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow

## This script should only be used on ICE

set -e # quits at first error

GOAL=$1
QUEMAX=$2
VERSION=genome #$3
OUT=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_${VERSION}
MODEL=$4
SYSTEM=$5 #smp, cluster, htc

RESULTS=${OUT}/results_sims_AJ_M${MODEL}

set -f

if [ -e switch${MODEL}.txt ] ; then

    IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/669a0e602306776ffa8d8be33e63574dfa2d1766/atmo_instance_ip.txt)

    echo ""
    echo "#################"
    date
    echo "pwd: ${PWD}"
    echo "Goal: ${GOAL}"
    echo "Que max: ${QUEMAX} "
    echo "Version: ${VERSION}"
    echo "Out path: ${OUT}"
    echo "Model: ${MODEL}"
    echo "System: ${SYSTEM}"

    if [ "$SYSTEM" == "ocelote" ] ; then
        qstat=/cm/shared/apps/pbspro/current/bin/qstat
        qsub=/cm/shared/apps/pbspro/current/bin/qsub
    else
        qstat=/usr/local/bin/qstat_local
        qsub=/usr/pbs/bin/qsub
    fi

    #check number of completed simulations
    echo "Check for ${GOAL} completed runs in $RESULTS"
#    COMP=$(ssh agladstein@${IP_ADDRESS} find /vol_c/results_macsSwig_AJmodels_${VERSION}/sim_values_AJ_M${MODEL} -type f | wc -l)
    COMP=5
    echo "${COMP} runs have completed"
    if [ "$COMP" -ge "$GOAL" ]; then
        echo "Goal completed"
        rm switch${MODEL}.txt
        echo "Goal completed. ${COMP} runs have completed in $RESULTS." | sendmail agladstein@email.arizona.edu
        exit
    else
        #check number of jobs in que
        if [ "$SYSTEM" == "cluster" ]; then
            JOBS=$($qstat | grep "agladstein" | grep -v smp | grep -v htc | cut -d " " -f1)
        elif [ "$SYSTEM" == "ocelote" ]; then
            JOBS=$($qstat | grep "agladstein" | cut -d " " -f1)
        else
            JOBS=$($qstat | grep "agladstein" | grep $SYSTEM | cut -d " " -f1)
        fi
        echo $JOBS
        n=0
        for j in $JOBS; do
	        q=$($qstat -t $j | grep -w "Q" | wc -l)
	        n=$(($n + $q))
        done
        echo "You have $n jobs in the que"
        if [ "$n" -ge "$QUEMAX" ]; then
	        echo "That's enough jobs in the que"
	        exit
        else
	        #create PBS scripts
            echo "./main_function_AJmodel_j2.sh ${SYSTEM} ${OUT} ${MODEL}"
            ./main_function_AJmodel_j2.sh ${SYSTEM} ${OUT} ${MODEL}

            cd /home/u15/agladstein/ABC
            echo "rsync -za macsSwig_AJmodels /xdisk/agladstein/macsSwig_AJmodels; cd /xdisk/agladstein/macsSwig_AJmodels"
            rsync -za macsSwig_AJmodels /xdisk/agladstein/
            cd /xdisk/agladstein/macsSwig_AJmodels

            echo "Submit to windfall"
            echo "$qsub HPC_workflow/PBS/macsargs_model${MODEL}_${SYSTEM}_windfall.pbs"
            $qsub HPC_workflow/PBS/macsargs_model${MODEL}_${SYSTEM}_windfall.pbs
        fi
    fi
else
    echo "switch${MODEL}.txt does not exist"
    exit
fi
