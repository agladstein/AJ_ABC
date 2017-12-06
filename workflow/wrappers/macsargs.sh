#!/bin/bash

set -e

MODEL=$1
JOB_ID=$2
SIM_SIZE=$3
NUM=$(echo $MODEL | cut -d "_" -f3)

# untar the model code
tar xzf model.tar.gz
cd model

# set up the venv
. /cvmfs/oasis.opensciencegrid.org/osg/sw/module-init.sh || true
module load python/2.7
# update the venv to the new location
virtualenv-2.7 workflow/macss_env
. workflow/macss_env/bin/activate

CMD="python gen_macsargs_${NUM}.py $JOB_ID $SIM_SIZE prior 0 ."
echo
echo "Running: $CMD"
$CMD

