#!/bin/bash                                                                                  

source /apps/root/6.18.04/bin/thisroot.csh 

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 20 20 6. ../../rootfiles/ca40_6gev_$i.root 50000 -P ../phase.txt
done

root -b -q '../GCF_Ca40_LUND.C(100,"ca40_6gev")'
