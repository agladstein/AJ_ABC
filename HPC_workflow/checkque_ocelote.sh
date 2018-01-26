#!/usr/bin/env bash
cd /home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow

## This script should only be used on Ocelote

set -e # quits at first error

GOAL=$1
QUEMAX=$2
VERSION=genome #$3
OUT=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_${VERSION}
MODEL=$4

RESULTS=${OUT}/results_sims_AJ_M${MODEL}

ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}

set -f

if [ -e switch${MODEL}.txt ] ; then

    IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/atmo_instance_ip.txt)

    echo ""
    echo "#################"
    date
    echo "pwd: ${PWD}"
    echo "Goal: ${GOAL}"
    echo "Que max: ${QUEMAX} "
    echo "Version: ${VERSION}"
    echo "Out path: ${OUT}"
    echo "Model: ${MODEL}"

    qstat=/cm/shared/apps/pbspro/current/bin/qstat
    qsub=/cm/shared/apps/pbspro/current/bin/qsub

    #check number of completed simulations
    echo "Check for ${GOAL} completed runs in $RESULTS"
    COMP_HPC=$(ssh agladstein@${IP_ADDRESS} find ${ATMO_DIR}/HPC/model${MODEL} -type f | wc -l)
    COMP_OSG=$(ssh agladstein@${IP_ADDRESS} find ${ATMO_DIR}/OSG -type f | xargs ssh agladstein@${IP_ADDRESS} cat | wc -l)
    COMP_CHTC=$(ssh agladstein@${IP_ADDRESS} find ${ATMO_DIR}/CHTC -type f | xargs ssh agladstein@${IP_ADDRESS} cat | wc -l)
    COMP=$(($COMP_CHTC + $COMP_HPC + $COMP_OSG))
    echo "${COMP} runs have completed"

    if [ "$COMP" -ge "$GOAL" ]; then
        echo "Goal completed"
        rm switch${MODEL}.txt
        echo "Goal completed. ${COMP} runs have completed in $RESULTS." | sendmail agladstein@email.arizona.edu
        exit
    else
        #check number of jobs in que
        JOBS=$($qstat | grep "agladstein" | cut -d " " -f1)
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
#            echo "./main_function_AJmodel_j2.sh ${OUT} ${MODEL}"
#            ./main_function_AJmodel_j2.sh ${OUT} ${MODEL}

            cd /home/u15/agladstein/ABC
            echo "rsync -za macsSwig_AJmodels /xdisk/agladstein/macsSwig_AJmodels; cd /xdisk/agladstein/macsSwig_AJmodels"
            rsync -za macsSwig_AJmodels /xdisk/agladstein/
            cd /xdisk/agladstein/macsSwig_AJmodels

            echo "Submit to standard"
            echo "$qsub HPC_workflow/PBS/macsargs_model${MODEL}_ocelote_standard.pbs"
            $qsub HPC_workflow/PBS/macsargs_model${MODEL}_ocelote_standard.pbs
        fi
    fi
else
    echo "switch${MODEL}.txt does not exist"
    exit
fi
