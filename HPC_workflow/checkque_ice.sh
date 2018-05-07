#!/usr/bin/env bash
cd /home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow

## This script should only be used on ICE

set -e # quits at first error

GOAL=$1
QUEMAX=$2
VERSION=update #$3
OUT=/rsgrps/mfh4/Ariella/macsSwig_AJmodels_${VERSION}
CHR=$3

RESULTS=${OUT}/results_sims_AJ_chr${CHR}

ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}/HPC/chr${CHR}

set -f

if [ -e switch${CHR}.txt ] ; then

    IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/atmo_instance_ip.txt)

    echo ""
    echo "#################"
    date
    echo "pwd: ${PWD}"
    echo "Goal: ${GOAL}"
    echo "Que max: ${QUEMAX} "
    echo "Version: ${VERSION}"
    echo "Out path: ${OUT}"
    echo "chr: ${CHR}"

    qstat=/cm/shared/apps/pbspro/current/bin/qstat
    qsub=/cm/shared/apps/pbspro/current/bin/qsub


    #check number of completed simulations
    echo "Check for ${GOAL} completed runs in $RESULTS"
    COMP=$(ssh agladstein@${IP_ADDRESS} find ${ATMO_DIR} -maxdepth 2 -type f | wc -l)
    echo "${COMP} runs have completed"
    if [ "$COMP" -ge "$GOAL" ]; then
        echo "Goal completed"
        rm switch${CHR}.txt
        echo "Goal completed. ${COMP} runs have completed in $RESULTS." | sendmail agladstein@email.arizona.edu
        exit
    else
        #check number of jobs in que
        JOBS=$($qstat | grep "agladstein" | cut -d " " -f1)
        echo $JOBS
        n=$(qstat -t -u agladstein | wc -l)
        echo "You have $n jobs"
        echo "$(qstat -t -u agladstein | grep -w "Q" | wc -l) are in the queue"
        echo "$(qstat -t -u agladstein | grep -w "R" | wc -l) are in running"
        if [ "$n" -ge "$QUEMAX" ]; then
	        echo "That's enough jobs"
	        exit
        else
	        #create PBS scripts
#            echo "./main_function_AJmodel_j2.sh ${SYSTEM} ${OUT} ${CHR}"
#            ./main_function_AJmodel_j2.sh ${SYSTEM} ${OUT} ${CHR}

            cd /home/u15/agladstein/ABC
            echo "rsync -za macsSwig_AJmodels /xdisk/agladstein/macsSwig_AJmodels; cd /xdisk/agladstein/macsSwig_AJmodels"
            rsync -za macsSwig_AJmodels /xdisk/agladstein/
            cd /xdisk/agladstein/macsSwig_AJmodels

            echo "Submit to windfall"
            echo "$qsub HPC_workflow/PBS/run_sims_update_chr${CHR}.pbs"
#            $qsub HPC_workflow/PBS/run_sims_update_chr${CHR}.pbs
        fi
    fi
else
    echo "switch${CHR}.txt does not exist"
    exit
fi
