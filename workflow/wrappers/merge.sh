#!/bin/bash

set -e

OUTFILE=$1 #output file
shift

# use the header from the first file
cp $1 $OUTFILE #the first input file becomes the top of the output file
shift
# for all the remaining, just the data
for FILE in "$@"; do
    tail -n +2 $FILE >>$OUTFILE
done



