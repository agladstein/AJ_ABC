GOAL=$1
RESULTS=$2
echo "Check for ${GOAL} completed runs in $RESULTS"
QUEMAX=$3
PBS=$4
qstat=/cm/shared/apps/pbspro/13.0.2.153173/bin/qstat
qsub=/cm/shared/apps/pbspro/13.0.2.153173/bin/qsub

#check number of completed simulations
COMP=$(ls $RESULTS/ | wc -l)
echo "${COMP} runs have completed"
if [ "$COMP" -ge "$GOAL" ]; then
    echo "Goal completed"
    crontab -r
else
    #check hrs left in group
    HRS=$(va | cut -f2 | tail -1 | cut -d ":" -f1)
    echo "${HRS} mfh hrs are left"
    DAYS=$(( $(($(cal | wc -w) - 9)) - $(date | cut -d " " -f3) ))
    echo "There are $DAYS days left in the month"
    BOUND=$(( $DAYS * 350))
    echo "You should leave $BOUND for the rest of the lab"
    if [ "$HRS" -le "$BOUND" ]; then
	echo "There are no hrs left to use"
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
	else
	    echo "$qsub $PBS"
	    $qsub $PBS
	fi
    fi
fi
    

