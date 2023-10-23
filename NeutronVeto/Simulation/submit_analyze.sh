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

#SBATCH --exclude=farm160122                                                           
#SBATCH --array=1-99

SEED=12345
TORUS=-1.0
#INPUT=/w/hallb-scshelf2102/clas/clase2/erins/repos/neutron-veto/Simulation_eN
PROGRAM_DIR=/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm_build/NeutronVeto
INPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/reconhipo_bkg
OUTPUT=/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/
NEVENTS=100
GCARD=rgm.gcard
YAML=rgm_mc.yaml


source /etc/profile.d/modules.sh
source /group/clas12/packages/setup.sh
module load clas12/pro #dev


# inputs: input hipo, output hipo, charge of particle, keep_good

##clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_neutrons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"good_neutrons_${SLURM_ARRAY_TASK_ID}\",0,1)" # look for good neutrons
#clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_neutrons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"bad_neutrons_${SLURM_ARRAY_TASK_ID}\",0,0)" # look for bad neutrons
#clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_neutrons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"fake_protons_${SLURM_ARRAY_TASK_ID}\",1,1)" # look for fake protons (2212 in e'n)


#clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_protons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"good_protons_${SLURM_ARRAY_TASK_ID}\",1,1)" # look for good protons
#clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_protons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"bad_protons_${SLURM_ARRAY_TASK_ID}\",1,0)" # look for bad protons
##clas12root -b -q "${INPUT}/N_skim.c(\"${INPUT}/eN_bknd/reconhipo_bkg/recon_protons_bkg_8_4_0_${SLURM_ARRAY_TASK_ID}.hipo\",\"fake_neutrons_${SLURM_ARRAY_TASK_ID}\",0,0)" # look for fake neutrons (2112 in e'p)


${PROGRAM_DIR}/N_getfeatures 0 ${OUTPUT}/ana_root/ana_neutrons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/ana_txt/ana_neutrons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/recon_neutrons_bkg_${SLURM_ARRAY_TASK_ID}.hipo



${PROGRAM_DIR}/N_getfeatures 1 ${OUTPUT}/ana_root/ana_protons_${SLURM_ARRAY_TASK_ID}.root ${OUTPUT}/ana_txt/ana_protons_${SLURM_ARRAY_TASK_ID}.txt ${INPUT}/recon_protons_bkg_${SLURM_ARRAY_TASK_ID}.hipo




# remember to check if charge is correct
