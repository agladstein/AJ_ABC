#!/usr/bin/env bash

set -e # quits at first error

cd /vol_c
tar cf - /vol_c | /vol_c/bin/drive push -exclude-ops delete,update -no-prompt -piped backup_atmo_vol_1T/backup_atmo_vol_1T_$(date +%m%d%Y%T).tar
