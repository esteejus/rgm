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


#SBATCH --array=5567,5568,5569,5570,5572,5573,5574,5575,5576,5577,5578,5579,5580,5581,5582,5583,5586,5587,5588,5589,5590,5591,5592,5593,5594,5595,5598,5599,5600,5601,5602,5603,5604,5606,5607,5608,5609,5610,5611,5612,5613,5614,5615,5616,5617,5618,5619,5620,5622,5623,5624,5625,5626,5627



INPUT=/lustre19/expphy/cache/clas12/rg-m/production/pass1/2gev/D/dst/recon
OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-efficiency/out_d_pfdn



source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev
module switch coatjava/8.4.0
module switch gemc/5.1
module load sqlite/dev



# Run d(e,e'pn) - RGM 2 GeV, p in FD
./neff_d_pfdn 0 2.07052 ${OUTPUT}/dFD_2gev_01${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/dFD_2gev_01${SLURM_ARRAY_TASK_ID}.pdf /w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronEfficiency/cutfile_d_neff_pfd.txt ${INPUT}/01${SLURM_ARRAY_TASK_ID}/rec_clas_01${SLURM_ARRAY_TASK_ID}.evio.*.hipo
