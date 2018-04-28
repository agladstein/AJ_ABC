#!/usr/bin/env bash

PBS=$1

RUN_DIR=/xdisk/agladstein/macsSwig_AJmodels
rsync -za ~/ABC/macsSwig_AJmodels /xdisk/agladstein/
cd ${RUN_DIR}

for i in {1..6}
do
    qsub ${PBS}
done