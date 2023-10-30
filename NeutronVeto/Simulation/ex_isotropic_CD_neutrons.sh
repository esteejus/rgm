#!/bin/bash                                                                                  

LUNDOUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/lundfiles
NEVENTS=10000

for i in {1..100}
do
    root -b -q "generate_neutrons.C(\"${LUNDOUT}/isotropic_neutrons_CD_$i.txt\",${NEVENTS},35,140)"
done

# should be 10000 events, 100 files
