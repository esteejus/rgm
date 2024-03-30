#!/bin/bash

# INDIR contains 100 e'p simulation files each with 10k events
# OUTDIR contains 100 background RG-A CLAS12 files each with 20k events

INPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/mchipo
BKGD_DIR=/cache/clas12/rg-m/production/bkgfiles/tor-1.00_sol-1.00/Cx4_5986MeV
OUTPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/mchipo_bkg

#rm ./ep_bkg/*

index=$1
nevents=10000

j=$(printf "%05d" $index)
bg-merger -d "BMT,BST,FMT,HTCC,DC,CND,CTOF,ECAL,FTOF" -n $nevents -b ${BKGD_DIR}/c_$j.hipo -i ${INPUT_DIR}/isotropic_protons_mc_${index}.hipo -o ${OUTPUT_DIR}/mc_protons_bkg_${index}.hipo
