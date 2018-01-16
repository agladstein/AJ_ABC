#!/usr/bin/env bash

set -e # quits at first error

VERSION=$1
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_${VERSION}

cd ${ATMO_DIR}
tar cf - ${ATMO_DIR} | /vol_c/bin/drive push -exclude-ops "delete,update" -no-prompt -piped backup_results_macsSwig_AJmodels_${VERSION}_$(date +%m%d%Y%T).tar
