#!/bin/bash                                                                                  

GCFGEN=/pathto/GCF_Generator_Suite/src/programs/genQE/
ROOTOUT=/root/output/dir/
LUNDOUT=/lund/output/dir/

for i in {1..1}
do
  $GCFGEN/genQE 6 6 5.98 $ROOTOUT/qe_d_6gev_$i.root 10 -P phase.txt
  root -b -q "GCF_to_LUND.C(\"${ROOTOUT}/qe_d_6gev_$i.root\",\"${LUNDOUT}/lund_qe_c_6gev_$i.txt\")"
done
