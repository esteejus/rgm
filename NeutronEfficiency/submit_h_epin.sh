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

#SBATCH --exclude=farm140235

#SBATCH --array=5019,5020,5022,5023,5024,5025,5026,5027,5028,5029,5030,5031,5032,5033,5034,5035,5036


INPUT=/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/H/dst/recon
OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-efficiency/out_h_epin
BUILD_DIR=/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm_build/NeutronEfficiency

source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
module switch coatjava/8.4.0
module switch gemc/5.1
module load sqlite/dev



# Run d(e,e'pn) - RGM 2 GeV, p in CD
${BUILD_DIR}/neff_h_epin 0 5.98636 ${OUTPUT}/hepin_6gev_01${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/hepin_6gev_01${SLURM_ARRAY_TASK_ID}.pdf /w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronEfficiency/cutfile_np_nosrc.txt ${INPUT}/01${SLURM_ARRAY_TASK_ID}/rec_clas_01${SLURM_ARRAY_TASK_ID}.evio.*.hipo
