#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=500
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_rgm                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=1:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err

#SBATCH --exclude=farm160122,farm180218,farm160138,farm140126,farm140236,farm160121,farm140215
                                                           
#SBATCH --array=1-99

SEED=12345
TORUS=-1.0
INPUT=/w/hallb-scshelf2102/clas/clase2/erins/repos/neutron-veto/Simulation_eN
OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/flatsim
NEVENTS=10000
GCARD=rgm.gcard
YAML=rgm_mc.yaml


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
module switch coatjava/8.4.0
module switch gemc/5.1
module load sqlite/dev


./N_getfeatures 0 ${OUTPUT}/sim_root/good_neutrons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/sim_txt/good_neutrons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/eN_bknd/reconhipo_bkg/recon_neutrons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo

./N_getfeatures 1 ${OUTPUT}/sim_root/fake_neutrons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/sim_txt/fake_neutrons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/eN_bknd/reconhipo_bkg/recon_protons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo


# remember to check if charge is correct
