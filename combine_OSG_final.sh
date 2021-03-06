#!/usr/bin/env bash

ATMO_PATH=$1

for dir in $ATMO_PATH/*
do
    workflow=${dir}/outputs
    if [ -e $workflow/sim_values.txt ] ; then
        echo $workflow
        MODEL=$(head -2 $workflow/sim_values.txt | tail -1 | cut -d '/' -f1 | cut -d '_' -f4)
        echo ${MODEL}

        cut -f2- $workflow/results_sims.txt > $workflow/cut_results_sims.txt

        if [ ${MODEL} == "M1" ] ; then
            NSTATS=$(head -1 $workflow/results_sims.txt | tr '\t' '\n' | wc -l)
            if [ -e input_ABC_OSG_${MODEL}_${NSTATS}.txt ] ; then
                paste $workflow/sim_values.txt $workflow/cut_results_sims.txt | grep -v "Asc_NAF" >>input_ABC_OSG_${MODEL}_${NSTATS}.txt
            else
                paste $workflow/sim_values.txt $workflow/cut_results_sims.txt >input_ABC_OSG_${MODEL}_${NSTATS}.txt
            fi

        else
            if [ -e input_ABC_OSG_${MODEL}.txt ] ; then
                paste $workflow/sim_values.txt $workflow/cut_results_sims.txt | grep -v "Asc_NAF" >>input_ABC_OSG_${MODEL}.txt
            else
                paste $workflow/sim_values.txt $workflow/cut_results_sims.txt >input_ABC_OSG_${MODEL}.txt
            fi
        fi
        rm $workflow/cut_results_sims.txt
    fi
done
