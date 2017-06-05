#!/usr/bin/env bash
cd /home/u15/agladstein/ABC/macsSwig_AJmodels

GOAL=$1
QUEMAX=$2
OUT=$3
MODEL=$4
SYSTEM=$5 #smp, cluster, htc

IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/8bfef8b8f3f7c43fd99832a323ef7130f98571bb/atmo_instance_ip.txt)

set -f

if [ -e switch${MODEL}.txt ] ; then

    echo "\n #################"
    date

    if [ "$SYSTEM" == "ocelote" ] ; then
        qstat=/cm/shared/apps/pbspro/current/bin/qstat
        qsub=/cm/shared/apps/pbspro/current/bin/qsub
    else
        qstat=/usr/local/bin/qstat_local
        qsub=/usr/pbs/bin/qsub
    fi

    RESULTS=${OUT}/results_sims_AJ_M${MODEL}
    date
    echo "Check for ${GOAL} completed runs in $RESULTS"

    #check number of completed simulations
    COMP=$(ssh agladstein@${IP_ADDRESS} find /vol_c/results_macsSwig_AJmodels_instant/sim_values_AJ_M${MODEL} -type f | wc -l)
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
        m=0
        p=0
        for j in $JOBS; do
	        q=$($qstat -t $j | grep -w "Q" | wc -l)
	        n=$(($n + $q))
	        r=$($qstat -t $j | grep "stan" | grep -w "R" | wc -l)
	        m=$(($m + $r))
	        s=$($qstat -t $j | grep "qual" | wc -l)
	        p=$(($p + $s))
        done
        echo "You have $n jobs in the que"
        if [ "$n" -ge "$QUEMAX" ]; then
	        echo "That's enough jobs in the que"
	        exit
        else
	        #create PBS scripts
            ./main_function_AJmodel_j2.sh ${SYSTEM} ${OUT} ${MODEL}

            #check standard hrs left in group
            SHRS=$(va | cut -f2 | tail -1 | cut -d ":" -f1)
	        DAYS=$(( $(($(cal | wc -w) - 9)) - $(($(date | cut -d " " -f3))) ))
	        SBOUND=$(( $DAYS * 1 + $n + $m))

	        echo "${SHRS} mfh standard hrs are left"
	        echo "There are $DAYS days left in the month"
            echo "You should leave $SBOUND for the rest of the lab"

            if [ "$SHRS" -le "$SBOUND" ]; then
                echo "There are no standard hrs left to use"

                if [ "$SYSTEM" == "smp" ] || [ "$SYSTEM" == "ocelote" ] ; then
                    #check qualified hrs left in group
                    QHRS=$(va | cut -f3 | tail -1 | cut -d ":" -f1)
                    QBOUND=$p
                    echo "${QHRS} mfh qualified hrs are left"
                    if [ "$QHRS" -gt "$QBOUND" ]; then
                        echo "Submit to qualified"
                        echo "$qsub model${MODEL}_${SYSTEM}_qualified.pbs"
                        $qsub model${MODEL}_${SYSTEM}_qualified.pbs
                        exit
                    else
                        echo "There are no qualified hrs left to use"
                    fi
                fi

                echo "Submit to windfall"
                echo "$qsub model${MODEL}_${SYSTEM}_windfall.pbs"
                $qsub model${MODEL}_${SYSTEM}_windfall.pbs

            else
                echo "Submit to standard"
	            echo "$qsub model${MODEL}_${SYSTEM}_standard.pbs"
	            $qsub model${MODEL}_${SYSTEM}_standard.pbs
	        fi
        fi
    fi
else
    exit
fi
