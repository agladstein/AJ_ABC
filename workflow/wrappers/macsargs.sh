#!/bin/bash

set -e

MODEL=$1
JOB_ID=$2
SIM_SIZE=$3
NUM=$(echo $MODEL | cut -d "_" -f3)

# untar the model code
tar xzf model.tar.gz
cd model

# set up the env (HOME is due to some CHTC machines not having our user)
export HOME=$PWD
export LD_LIBRARY_PATH=$PWD/workflow/macss_env/lib
export PATH=$PWD/workflow/macss_env/bin:$PATH
which -a python

CMD="python gen_macsargs_${NUM}.py $JOB_ID $SIM_SIZE prior 0 ."
echo
echo "Running: $CMD"
$CMD
mv macsargs_${JOB_ID}.txt ../
