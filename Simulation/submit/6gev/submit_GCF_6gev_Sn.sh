#!/bin/bash                                                                                  

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 50 70 6. ../../rootfiles/sn_6gev_$i.root 50000 -P ../phase.txt
done

root -b -q '../GCF_multiFoil_Sn_LUND.C(100,"sn_6gev")'
