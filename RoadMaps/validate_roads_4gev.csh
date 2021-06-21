#!/bin/tcsh -f

set FILE=$1
set VALIDATE_FILE=/volatile/clas12/rg-b/production/recon/pass0/florian/calib/recon/011299/rec_clas_011299.evio.00715-00719.hipo
#set VALIDATE_FILE=/lustre19/expphy/cache/clas12/rg-b/production/recon/pass0/v24.3/dst/train/inc/inc_011299.hipo
#set VALIDATE_FILE=/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/recon/006156/rec_clas_006156.evio.00000-00004.hipo
source /group/clas12/packages/setup.csh
module load clas12/pro

/work/clas12/users/devita/roads/coatjava-roads/bin/dict-validate -dict ${FILE} -i ${VALIDATE_FILE} -wire 1 -mode 0 -dupli 1 -charge -1 -pid 11 -n 100000 -sector 0
