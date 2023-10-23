#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=recon_bkg
                                                                 
#SBATCH --partition=production                                                               
#SBATCH --time=72:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-100                                                                                                     

SEED=12345
TORUS=-1.0

INPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/mchipo_bkg
OUTPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/reconhipo_bkg
SUPPORT_DIR=/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronVeto/Simulation

OUTPUT_FILE=neutrons_bkg
#OUTPUT_FILE=protons_bkg
NEVENTS=10000
GCARD=rgm.gcard
YAML=rgm_mc.yaml


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
#module switch coatjava/8.4.0
module switch coatjava/9.0.0
#module switch gemc/5.1 # use 5.3
module load sqlite/dev


recon-util -n $NEVENTS -y ${SUPPORT_DIR}/$YAML -i ${INPUT_DIR}/mc_${OUTPUT_FILE}_${SLURM_ARRAY_TASK_ID}.hipo -o ${OUTPUT_DIR}/recon_${OUTPUT_FILE}_${SLURM_ARRAY_TASK_ID}.hipo
