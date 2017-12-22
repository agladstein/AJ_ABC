#!/bin/bash

set -e

MODEL=$1
JOB_ID=$2
INPUT_FILE=$3
SIM_SIZE=$4

# untar the model code
tar xzf model.tar.gz
cd model

# set up the env (HOME is due to some CHTC machines not having our user)
export HOME=$PWD
export LD_LIBRARY_PATH=$PWD/workflow/macss_env/lib
export PATH=$PWD/workflow/macss_env/bin:$PATH
which -a python

CMD="python $MODEL $JOB_ID $INPUT_FILE $SIM_SIZE 0 prior 0"
echo
echo "Running: $CMD"
$CMD
mv results* ../
mv *.match ../


