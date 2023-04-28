#!/bin/bash                                                                                  

source /site/12gev_phys/softenv.sh 2.4
source /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/geant4/4.10.06.p02MT/bin/geant4.sh

GCFGEN=../../GCF_build/programs/genQE
ROOTOUT=rootfiles
LUNDOUT=lundfiles

for i in {1..100}
do
  $GCFGEN/genQE 6 6 6.0 $ROOTOUT/qe_c_6gev_$i.root 10000 -P phase.txt
  root -b -q "GCF_to_LUND.C(\"${ROOTOUT}/qe_c_6gev_$i.root\",\"${LUNDOUT}/lund_qe_c_6gev_$i.txt\")"
done
