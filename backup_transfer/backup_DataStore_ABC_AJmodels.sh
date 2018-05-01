#!/usr/bin/env bash

set -e # quits at first error

VERSION=$1

# copy current directory to Data Store
ATMO_PATH=/vol_c/ABC_AJmodels_${VERSION}
IPLANT_PATH=/iplant/home/agladstein/ABC_AJmodels

echo 'iput -K -f -r ${ATMO_PATH} ${IPLANT_PATH}'
iput -K -f -r ${ATMO_PATH} ${IPLANT_PATH}

### backup current directory as gzip if directory changed
#ibun ?