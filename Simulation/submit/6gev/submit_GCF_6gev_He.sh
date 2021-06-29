#!/bin/bash                                                                                  

source /apps/root/6.18.04/bin/thisroot.csh 

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 2 2 6. ../../rootfiles/he_6gev_$i.root 50000 -P ../phase.txt
done

root -b '../GCF_He_LUND.C(100,"he_6gev")'
