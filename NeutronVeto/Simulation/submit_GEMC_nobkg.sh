#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=gemc_nobkg                                                                                            
#SBATCH --partition=production                                                               
#SBATCH --time=72:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
#SBATCH --array=1-100                                                                                                     

SEED=12345
TORUS=-1.0

INPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/lundfiles
OUTPUT_DIR=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/mchipo
SUPPORT_DIR=/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronVeto/Simulation

#OUTPUT_FILE=neutrons_bkg
#OUTPUT_FILE=protons_bkg

#OUTPUT_FILE=isotropic_neutrons
OUTPUT_FILE=isotropic_protons

NEVENTS=10000
GCARD=rgm.gcard
YAML=rgm_mc.yaml


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
#module switch coatjava/8.4.0
#module switch coatjava/8.6.0 # want to use current version (10.0.2?) 
module switch gemc/5.2 # switched to 5.3 from 5.1 (has latest cad drawings)
module load sqlite/dev


gemc -USE_GUI=0  -RANDOM=${SEED} -SCALE_FIELD="TorusSymmetric, $TORUS" -SCALE_FIELD="clas12-newSolenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ${INPUT_DIR}/${OUTPUT_FILE}_CD_${SLURM_ARRAY_TASK_ID}.txt" -OUTPUT="hipo, ${OUTPUT_DIR}/${OUTPUT_FILE}_mc_${SLURM_ARRAY_TASK_ID}.hipo" ${SUPPORT_DIR}/$GCARD
