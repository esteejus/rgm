#!/bin/bash                                                                                  

source /apps/root/6.18.04/bin/thisroot.csh 

for i in {1..100}
do
../../GCF_Generator_Suite/programs/genQE/genQE 1 1 2. ../../rootfiles/d_2gev_$i.root 50000 -P ../phase.txt
done

root -b '../GCF_D_LUND.C(100,"d_2gev")'
