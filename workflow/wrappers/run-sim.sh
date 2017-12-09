#!/bin/bash

set -e

MODEL=$1
CHR=$2
MACS_ARGS_FILE=$3
SNP_FILE=$4

# untar the model code
tar xzf model.tar.gz
cd model

# set up the env
export LD_LIBRARY_PATH=$PWD/workflow/macss_env/lib
export PATH=$PWD/workflow/macss_env/bin:$PATH
which -a python

mv ../macsargs_*.txt .

CMD="python $MODEL $CHR $MACS_ARGS_FILE $SNP_FILE 0 ."
echo
echo "Running: $CMD"
$CMD
mv results* ../
mv *.match ../


