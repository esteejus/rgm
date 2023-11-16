#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=500
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=veto_sim
#SBATCH --partition=production                                                               
#SBATCH --time=3:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err

#SBATCH --exclude=farm160122
                                                           
#SBATCH --array=1-99

SEED=12345
TORUS=-1.0
INPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/reconhipo_bkg
OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev


./N_getfeatures 0 ${OUTPUT}/ana_root/neutrons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/ana_txt/neutrons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/recon_neutrons_bkg_${SLURM_ARRAY_TASK_ID}.hipo

./N_getfeatures 1 ${OUTPUT}/ana_root/protons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/ana_txt/protons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/recon_protons_bkg_${SLURM_ARRAY_TASK_ID}.hipo


# remember to check if charge is correct
