#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_genie_c12_4gev                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=30:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-100                                                                                                     

NEVENTS=10000
TORUS=-1.0
FILE_PREFIX=c12_4gev #Change file prefix for your simulation

GCARD=./rgm.gcard
YAML=./rgm_mc.yaml

source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load sqlite/dev
module load clas12/dev

gemc -USE_GUI=0  -SCALE_FIELD="TorusSymmetric, $TORUS" -SCALE_FIELD="clas12-newSolenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ../lundfiles/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.dat" -OUTPUT="hipo, ../mchipo/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD
recon-util -y $YAML -i ../mchipo/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo -o ../reconhipo/recon_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo
