#!/usr/bin/env bash

set -oe

DATA_DIR=$1
CHR=$2

CGI=ALL.atDNA.biAllelicSNPnoDI.genotypes_maf0.005_YRI9_CEU9_CHB4_chr${CHR}.tped
CGIarray=ALL.atDNA.biAllelicSNPnoDI.genotypes_maf0.005_YRI9_CEU9_CHB4_hg18_Behar_HGDP_FtDNA_chr${CHR}.tped
ARRAY=Behar_HGDP_FtDNA_Jews_MidEast_subset_chr${CHR}.tped
GERMLINE=Behar_HGDP_FtDNA_Jews_MidEast_subset_YRI9_CEU9_CHB4_chr${CHR}.ped

/vol_c/env/simprily_env2.7/bin/python main_function_realdata_M23.py ${CHR} ${DATA_DIR} ${CGI} ${CGIarray} ${ARRAY} ${GERMLINE}