#!/bin/bash

set -e

MODEL=$1
JOB_ID=$2
INPUT_FILE=$3
SIM_SIZE=$4

# untar the model code
tar xzf model.tar.gz
cd model

# set up the venv
. /cvmfs/oasis.opensciencegrid.org/osg/sw/module-init.sh || true
module load python/2.7
# update the venv to the new location
virtualenv-2.7 workflow/macss_env
. workflow/macss_env/bin/activate

CMD="python $MODEL $JOB_ID $INPUT_FILE $SIM_SIZE 0 prior 0"
echo
echo "Running: $CMD"
$CMD

# output files need unique names in the top level dir
tar czf ../outputs-$JOB_ID.tar.gz results_sims_* sim_values_*

