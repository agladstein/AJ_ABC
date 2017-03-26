#!/usr/bin/env bash
cd /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels

GOAL=$1
RESULTS=$2
echo "Check for ${GOAL} completed runs in $RESULTS"
QUEMAX=$3
PBS=$4
SYSTEM=$5 #smp, clu, htc
qstat=/usr/local/bin/qstat_local
qsub=/usr/pbs/bin/qsub

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
	    #check hrs left in group
	    QHRS=$(va | cut -f3 | tail -1 | cut -d ":" -f1)
	    QBOUND=30000 #0
	    SHRS=$(va | cut -f2 | tail -1 | cut -d ":" -f1)
	    DAYS=$(( $(($(cal | wc -w) - 9)) - $(date | cut -d " " -f3) ))
	    SBOUND=$(( $DAYS * 350))

	    echo "${QHRS} mfh qualified hrs are left"
	    if [ "$QHRS" -gt "$QBOUND" ]; then
	        echo "Submit to qualified"
	    else
	        echo "There are no qualified hrs left to use"
	        echo "Try standard"


#	    echo "${HRS} mfh standard hrs are left"
#	    echo "There are $DAYS days left in the month"
#	    echo "You should leave $BOUND for the rest of the lab"
#
#	if [ "$SHRS" -le "$SBOUND" ]; then
#	    echo "There are no hrs left to use"
#	    echo "Submit to windfall"
#	else
#
#	    echo "Submit to standard"
#	    echo "$qsub $PBS"
#	    #$qsub $PBS
#
	    fi
    fi
fi
    

