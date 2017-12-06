#!/bin/bash

set -e

MODEL=$1
CHR=$2
MACS_ARGS_FILE=$3
SNP_FILE=$4

# untar the model code
tar xzf model.tar.gz
cd model

# set up the venv
. /cvmfs/oasis.opensciencegrid.org/osg/sw/module-init.sh || true
module load python/2.7
# update the venv to the new location
virtualenv-2.7 workflow/macss_env
. workflow/macss_env/bin/activate

CMD="python $MODEL $CHR $MACS_ARGS_FILE $SNP_FILE 0 ."
echo
echo "Running: $CMD"
$CMD



