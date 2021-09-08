#!/bin/bash                                                                                  

source /apps/root/6.18.04/bin/thisroot.csh 

for i in {1..100}
do
../../GCF_Generator_Suite/programs/genQE/genQE 6 6 4. ../../rootfiles/c_4gev_$i.root 50000 -P ../phase.txt
done

root -b '../GCF_singleFoil_C_LUND.C(100,"c_4gev")'
