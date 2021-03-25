#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=test                                                                                              
#SBATCH --partition=production                                                                                        
#SBATCH --mail-user=esteejus@mit.edu                                                                                 
#SBATCH --time=4:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-50                                                                                                     

source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro

cd /work/clas12/users/esteejus/clas12analysis/ForFarm/LowEnergyDeuteron

srun gemc -USE_GUI=0 -N=20000 -INPUT_GEN_FILE="lund, ./lundfiles/lund_qe_$SLURM_ARRAY_TASK_ID.dat" -OUTPUT="evio, ./eviofiles/out_$SLURM_ARRAY_TASK_ID.ev" /group/clas12/gemc/4.4.1/config/rgb_spring2019.gcard
srun evio2hipo -t -1 -s -1 -i ./eviofiles/out_$SLURM_ARRAY_TASK_ID.ev -o ./mchipo/qe_$SLURM_ARRAY_TASK_ID.hipo
srun recon-util -y /group/clas12/gemc/4.4.1/config/rgb_spring2019.yaml -i ./mchipo/qe_$SLURM_ARRAY_TASK_ID.hipo -o ./reconhipo/recon_qe_$SLURM_ARRAY_TASK_ID.hipo
