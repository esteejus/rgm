#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_rgm                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=72:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-100                                                                                                     


TORUS=-1.0
RUN=flatgen
NEVENTS=10000
GCARD=submit/rgm.gcard
YAML=submit/rgm_mc.yaml

source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/dev
module switch gemc/5.1
module load sqlite/dev

time gemc -USE_GUI=0  -SCALE_FIELD="TorusSymmetric, $TORUS" -SCALE_FIELD="clas12-newSolenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, lundfiles/${RUN}_${SLURM_ARRAY_TASK_ID}.txt" -OUTPUT="hipo, mchipo/${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD
time recon-util -n $NEVENTS -y $YAML -i mchipo/${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo -o reconhipo/recon_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo

