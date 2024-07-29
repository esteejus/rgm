#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --account=clas12
#SBATCH --job-name=mc_rgm_gcf
#SBATCH --partition=production
#SBATCH --time=10:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --array=1-401 #Number of files 1-N

#ADD YOUR PATHS FROM INPUT TO OUTPUT HERE
#exacutables
GCF_EX=/w/hallb-scshelf2102/clas12/users/awild/GCF_Generators/GCF_Generator_Suite/src/programs/genQE/genQE
LUND_EX=/w/hallb-scshelf2102/clas12/users/awild/RGM/rgm/Simulation/GCF_to_LUND.C

#output
#change this to your specific output path
OUTPATH=/volatile/clas12/rg-m/awild/mc
FILE_PREFIX=src_qe_he_6gev

#input
#Here is where I put the configurations files I want for this simulation
INPATH=/w/hallb-scshelf2102/clas12/users/awild/RGM/CLAS_note/Simulation
PHASE_SPACE=${INPATH}/phase.txt
GCARD=${INPATH}/rgm_fall2021_He.gcard
YAML=${INPATH}/rgm_fall2021-ai_6Gev.yaml

Z=2
N=2
BEAM_E=5.98636  #5.98636, 4.02962, 2.07052
NEVENTS=20000  #20000
#-1.0 for inbending(6,4 GeV) 0.5 for outbending (2 Gev)
TORUS=-1.0 
TARGET=liquid #Targets: liquid, 4-foil, 1-foil, Ar, Ca
SIGMACM=0.100 #GeV/c


#DON'T NEED TO TOUCH BELOW HERE UNLESS YOU NEED TO
ROOTOUT=$OUTPATH/rootfiles
LUNDOUT=${OUTPATH}/lundfiles
MCOUT=${OUTPATH}/mchipo
RECONOUT=${OUTPATH}/reconhipo

#GCF Generator
$GCF_EX $Z $N $BEAM_E $ROOTOUT/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root $NEVENTS -s $SIGMACM -P $PHASE_SPACE -v -C -O

#TO LUND File
root -b -q "${LUND_EX}(\"${ROOTOUT}/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root\",\"${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt\",\"${TARGET}\")"

#SUBMIT GEMC MC
gemc -USE_GUI=0  -SCALE_FIELD="binary_torus, $TORUS" -SCALE_FIELD="binary_solenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt" -OUTPUT="hipo, ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD

#RECONSTRUCTION
recon-util -y $YAML -n $NEVENTS -i ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo -o ${RECONOUT}/recon_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo
