#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=500                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_rgm_gcf                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=00:30:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-1000 #Number of files 1-N                                                                                                

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
GCFGEN=/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Simulation/GCF_Generator_Suite/src/programs/genQE/ #path to compiled GCF generator
OUTPATH=/volatile/clas12/rg-m/mc/

Z=2
N=2
BEAM_E=5.98636  #5.98636, 4.02962, 2.07052
NEVENTS=10000
FILE_PREFIX=qe_he_6gev_sigmacm_130
TARGET=liquid #Targets: liquid, 4-foil, 1-foil, Ar, Ca
SIGMACM=0.130 #GeV/c

#DON'T NEED TO TOUCH BELOW HERE UNLESS YOU NEED TO
ROOTOUT=$OUTPATH/rootfiles/  #GCF root file output path                        
LUNDOUT=$OUTPATH/lundfiles/  #LUND file output path

#GCF Generator

$GCFGEN/genQE $Z $N $BEAM_E $ROOTOUT/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root $NEVENTS -s $SIGMACM -P phase.txt
#TO LUND File
root -b -q "GCF_to_LUND.C(\"${ROOTOUT}/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root\",\"${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt\",\"${TARGET}\")"
