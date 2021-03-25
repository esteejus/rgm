#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=512
#SBATCH --account=clas12
#SBATCH --job-name=test
#SBATCH --partition=production
#SBATCH --mail-user=esteejus@mit.edu
#SBATCH --time=2:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=//farm_out/%u/%x-%j-%N.err
#SBATCH --array=0-100

source /etc/profile.d/modules.sh   
source /group/clas12/packages/setup.sh
module load clas12/pro
module switch clas12root/1.5.1

cd /work/clas12/users/esteejus/clas12analysis/ForFarm/LowEnergyDeuteron

FILES=(/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/inc/*)

srun -W 20 clas12root -l -b HighEnergyReader.C+ --in=${FILES[$SLURM_ARRAY_TASK_ID]}
