#!/bin/bash                                                                                   

for i in {1..100}
do
  ../../GCF_Generator_Suite/programs/genQE/genQE 20 28 6. ../../rootfiles/ca48_6gev_$i.root 50000 -P ../phase.txt
done

root -b -q '../GCF_Ca48_LUND.C(100,"ca48_6gev")'
