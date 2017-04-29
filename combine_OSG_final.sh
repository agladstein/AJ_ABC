#!/usr/bin/env bash

for workflow in macsswig_simsaj_*
do
    if [ -e $workflow/sim_values.txt ] ; then
        echo $workflow
        MODEL=$(head -2 $workflow/sim_values.txt | tail -1 | cut -d '/' -f1 | cut -d '_' -f4)
        echo ${MODEL}
        if [ -e sim_values_OSG_${MODEL}.txt ] ; then
            echo 'grep -v "Asc_NAF" $workflow/sim_values.txt >>sim_values_OSG_${MODEL}.txt'
            grep -v "Asc_NAF" $workflow/sim_values.txt >>sim_values_OSG_${MODEL}.txt
        else
            echo 'cat $workflow/sim_values.txt >sim_values_OSG_${MODEL}.txt'
            cat $workflow/sim_values.txt >sim_values_OSG_${MODEL}.txt
        fi
        if [ -e results_sims_OSG_${MODEL}.txt ] ; then
            echo 'grep -v "SegS_Af_CGI" $workflow/results_sims.txt >>results_sims_OSG_${MODEL}.txt'
            grep -v "SegS_Af_CGI" $workflow/results_sims.txt >>results_sims_OSG_${MODEL}.txt
        else
            echo 'cat $workflow/results_sims.txt >results_sims_OSG_${MODEL}.txt'
            cat $workflow/results_sims.txt >results_sims_OSG_${MODEL}.txt
        fi
    fi
done