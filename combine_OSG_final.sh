#!/usr/bin/env bash

for workflow in macsswig_simsaj_*
do
    if [ -e $workflow/sim_values.txt ] ; then
        echo $workflow
        MODEL=$(head -2 $workflow/sim_values.txt | tail -1 | cut -d '/' -f1 | cut -d '_' -f4)
        echo ${MODEL}

        if [ -e input_ABC_OSG_${MODEL}.txt ] ; then
            paste $workflow/sim_values.txt $workflow/results_sims.txt | grep -v "Asc_NAF" >>input_ABC_OSG_${MODEL}.txt
        else
            paste $workflow/sim_values.txt $workflow/results_sims.txt >input_ABC_OSG_${MODEL}.txt
        fi
    fi
done