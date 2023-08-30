#!/bin/bash                                                                                  

#Code runs GCF genrator
#Generates a root file ouput from GCF and saves in ROOTOUT path
#Runs GCF_to_LUND.C macro which converts the root files into LUND (.txt) files for input into GEMC
#Specify to GCF_to_LUND.C which target you want to randomly generate over their extents
#liquid - standard liquid cell
#4-foil - 4-foil soild target flags
#1-foil - single foil used +2.5cm upstream
#Ar - LArgon cell
#Ca - Calcium targets

#Target options: liquid, 4-foil, 1-foil, Ar, Ca

source environment_gemc.sh

#ADD YOUR PATHS FROM INPUT TO OUTPUT HERE
GCFGEN=/pathto/GCF_Generator_Suite/src/programs/genQE/ #path to compiled GCF generator
ROOTOUT=/root/output/dir/  #GCF root file output path                        
LUNDOUT=/lund/output/dir/  #LUND file output path

Z=6
N=6
BEAM_E=5.98636  #5.98636, 4.02962, 2.07052
NEVENTS=10
NFILES=1
FILE_PREFIX=qe_d_6gev
TARGET=liquid #Targets: liquid, 4-foil, 1-foil, Ar, Ca


#DON'T NEED TO TOUCH BELOW HERE UNLESS YOU NEED TO
for (( i = 1; i <= $NFILES; i++ )) 
do
  $GCFGEN/genQE $Z $N $BEAM_E $ROOTOUT/${FILE_PREFIX}_$i.root $NEVENTS -P phase.txt
  root -b -q "GCF_to_LUND.C(\"${ROOTOUT}/${FILE_PREFIX}_$i.root\",\"${LUNDOUT}/lund_${FILE_PREFIX}_$i.txt\",\"${TARGET}"\)"
done
