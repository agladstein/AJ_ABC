#!/usr/bin/env bash

VERSION=$1
OSG_DIR=/local-scratch/agladstein/workflows/macsswig_simsaj_$VERSION
ATMO_DIR=/vol_c/results_macsSwig_AJmodels_$VERSION
echo ${ATMO_DIR}

IP_ADDRESS=$(curl https://gist.githubusercontent.com/agladstein/2bdc122f50314f2a4c7cbc9544e7a325/raw/8bfef8b8f3f7c43fd99832a323ef7130f98571bb/atmo_instance_ip.txt)

echo ssh agladstein@${IP_ADDRESS} mkdir -p $ATMO_DIR/OSG
ssh agladstein@${IP_ADDRESS} mkdir -p $ATMO_DIR/OSG

for w in $OSG_DIR*
do
    if [ -d $w/outputs ]; then
	BUCKET=$(echo $w | cut -d "_" -f4)
        echo $w/outputs
	echo rsync -avz $w/outputs agladstein@${IP_ADDRESS}:$ATMO_DIR/OSG/${BUCKET}
	rsync -avz $w/outputs agladstein@${IP_ADDRESS}:$ATMO_DIR/OSG/${BUCKET}
    fi
done




