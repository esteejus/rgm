#!/bin/bash                                                                                  

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 6 6 6. ../../rootfiles/c_6gev_$i.root 50000 -P ../phase.txt
done

root -b -q '../GCF_multiFoil_C_LUND.C(1000,"c_6gev")'
