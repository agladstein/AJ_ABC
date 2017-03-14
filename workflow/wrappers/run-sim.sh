#!/bin/bash

set -e

JOB_ID=$1
INPUT_FILE=$2
SIM_SIZE=$3

# untar the model code
tar xzf model.tar.gz
cd model

# set up the venv
. /cvmfs/oasis.opensciencegrid.org/osg/sw/module-init.sh || true
module load python/2.7
# update the venv to the new location
virtualenv-2.7 workflow/macss_env
. workflow/macss_env/bin/activate

CMD="python run_sims_AJmodel1_chr1_all.py $JOB_ID $INPUT_FILE $SIM_SIZE 0 prior 0"
echo
echo "Running: $CMD"
$CMD

# output files need unique names in the top level dir
tar czf ../outputs-$JOB_ID.tar.gz results_sims_AJ_M1 sim_values_AJ_M1

