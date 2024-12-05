#!/bin/bash                                                                                                          
#SBATCH --nodes=1                                                                                                    
#SBATCH --ntasks=1                                                                                                   
#SBATCH --mem-per-cpu=700
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=veto_goodn
#SBATCH --partition=production                                
#SBATCH --time=3:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err


#SBATCH --array=5567,5568,5569,5570,5572,5573,5574,5575,5576,5577,5578,5579,5580,5581,5583,5586,5587,5588,5589,5590,5591,5592,5593,5594,5595,5598,5599,5600,5601,5602,5603,5604,5606,5607,5608,5609,5610,5611,5612,5613,5614,5615,5616,5617,5618,5619,5620,5622,5623,5624,5625,5626,5627



OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/d_2gev
BUILD_DIR=/w/hallb-scshelf2102/clas12/erins/build_rgm/NeutronVeto


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
#module load clas12/pro #dev
module load sqlite/dev

# is energy = not-int a problem???


# COOKED RG-M

INPUT=/lustre19/expphy/cache/clas12/rg-m/production/pass1/2gev/D/dst/recon/

${BUILD_DIR}/D_getfeatures 2.07052 1 ${OUTPUT}/2gev_root/goodn_e5_pCD_01${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/2gev_txt/goodn_e5_pCD_01${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/01${SLURM_ARRAY_TASK_ID}/rec_clas_01${SLURM_ARRAY_TASK_ID}.evio.*.hipo # look for good neutrons

