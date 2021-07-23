#!/bin/bash                                                                                   

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 18 22 6. ../../rootfiles/ar_6gev_$i.root 50000 -P ../phase.txt
done

root -b -q '../GCF_Ar_LUND.C(100,"ar_6gev")'
