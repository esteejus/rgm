#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_rgm                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=20:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-100                                                                                                     

TORUS=0.5
RUN=d_2gev

source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro
module switch clas12root/1.7.2

srun gemc -USE_GUI=0  -SCALE_FIELD="TorusSymmetric, $TORUS" -SCALE_FIELD="clas12-newSolenoid, -1.0" -N=50000 -INPUT_GEN_FILE="lund, ../../lundfiles/lund_qe_${RUN}_${SLURM_ARRAY_TASK_ID}.dat" -OUTPUT="evio, ../../eviofiles/out_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.ev" /group/clas12/gemc/4.4.1/config/rgb_spring2019.gcard
srun evio2hipo -t $TORUS -s -1 -i ../../eviofiles/out_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.ev -o ../../mchipo/qe_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo
srun recon-util -y /group/clas12/gemc/4.4.1/config/rgb_spring2019.yaml -i ../../mchipo/qe_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo -o ../../reconhipo/recon_qe_${RUN}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo
