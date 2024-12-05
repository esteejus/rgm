#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=700
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=veto_ppipn
#SBATCH --partition=production                                
#SBATCH --time=3:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err


                                          
#SBATCH --array=5045,5046,5047,5049,5050,5051,5052,5053,5054,5055,5056,5057,5058,5059,5060,5061,5062,5065,5066,5067,5072,5073,5074,5075,5077,5078,5079,5081,5082,5093,5094,5095,5096,5097,5098,5099,5100,5101,5102,5103,5104,5105,5106,5434,5435,5436,5437,5439,5441,5442,5443,5444,5445,5447,5448,5449,5450,5451,5452,5454,5455,5456



OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/d_6gev
BUILD_DIR=/w/hallb-scshelf2102/clas12/erins/build_rgm/NeutronVeto



source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
#module switch coatjava/8.4.0
module load sqlite/dev

# is energy = not-int a problem???


# COOKED RG-M

INPUT=/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon

${BUILD_DIR}/D_getfeatures_ppim 5.98636 0 ${OUTPUT}/6gev_root/ppipn_e5_CD_01${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/6gev_txt/ppipn_e5_CD_01${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/01${SLURM_ARRAY_TASK_ID}/rec_clas_01${SLURM_ARRAY_TASK_ID}.evio.*.hipo # look for bad neutrons using d(e,e'p pi-)p

