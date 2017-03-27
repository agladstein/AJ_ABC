#!/usr/bin/env bash
cd /home/u15/agladstein/ABC/macsSwig_AJmodels

GOAL=$1
QUEMAX=$2
OUT=$3
MODEL=$4
SYSTEM=$5 #smp, clu, htc
qstat=/usr/local/bin/qstat_local
qsub=/usr/pbs/bin/qsub

RESULTS=${OUT}/results_sims_AJ_M${MODEL}
echo "Check for ${GOAL} completed runs in $RESULTS"

#check number of completed simulations
COMP=$(ls $RESULTS/ | wc -l)
echo "${COMP} runs have completed"
if [ "$COMP" -ge "$GOAL" ]; then
    echo "Goal completed"
    #crontab -r
else
    #check number of jobs in que	
    JOBS=$($qstat | grep "agladstein" | grep $SYSTEM | cut -d " " -f1)
    echo $JOBS
    n=0
    for j in $JOBS; do
	q=$($qstat -t $j | grep -w "Q" | wc -l)
	n=$(($n + $q))
    done
    echo "You have $n jobs in the que"
    if [ "$n" -ge "$QUEMAX" ]; then
	    echo "That's enough jobs in the que"

    else
	    #create PBS scripts
	    if [ "$SYSTEM" == "clu" ]; then
            ./main_function_AJmodel_j2.sh cluster /rsgrps/mfh4/Ariella/macsSwig_AJmodels ${MODEL}
        else
            ./main_function_AJmodel_j2.sh ${SYSTEM} /rsgrps/mfh4/Ariella/macsSwig_AJmodels ${MODEL}
        fi

	    #check hrs left in group
	    QHRS=$(va | cut -f3 | tail -1 | cut -d ":" -f1)
	    QBOUND=0
	    SHRS=$(va | cut -f2 | tail -1 | cut -d ":" -f1)
	    DAYS=$(( $(($(cal | wc -w) - 9)) - $(date | cut -d " " -f3) ))
	    SBOUND=$(( $DAYS * 350))

	    echo "${QHRS} mfh qualified hrs are left"
	    if [ "$QHRS" -gt "$QBOUND" ]; then
	        echo "Submit to qualified"
	        echo "qsub model${MODEL}_${SYSTEM}_qualified.pbs"
	    else
	        echo "There are no qualified hrs left to use"

	        echo "Try standard"
	        echo "${SHRS} mfh standard hrs are left"
	        echo "There are $DAYS days left in the month"
	        echo "You should leave $SBOUND for the rest of the lab"

        	if [ "$SHRS" -le "$SBOUND" ]; then
	            echo "There are no standard hrs left to use"
	            echo "Submit to windfall"
	        else
    	    echo "Submit to standard"
#	        echo "$qsub $PBS"
#	        $qsub $PBS
            fi
	    fi
    fi
fi
    

