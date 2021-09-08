#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1500
#SBATCH --account=clas12
#SBATCH --job-name=rgm_roads_4_4gev_fastmc
#SBATCH --partition=production
#SBATCH --time=20:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=//farm_out/%u/%x-%j-%N.err
#SBATCH --array=0-100


source /group/clas12/packages/setup.sh
module load clas12/pro


RANDOM=$(date +%s%N | cut -b10-19)
SEED=$(( $RANDOM % 3000 + 1 ))
TORUS=-1.0
SOLENOID=-1.0
CHARGE=-1
EVENTS=1000000

MIN_MOM=0.3
MAX_MOM=4.4
MIN_PHI=-60
MAX_PHI=60
MIN_THETA=5
MAX_THETA=40
MIN_VZ=-10
MAX_VZ=10

/work/clas12/users/devita/roads/coatjava-roads/bin/dict-maker -t $TORUS -s $SOLENOID -q $CHARGE -n $EVENTS -thmin $MIN_THETA -thmax $MAX_THETA -pmin $MIN_MOM -pmax $MAX_MOM -phimin $MIN_PHI -phimax $MAX_PHI -seed $((SEED)) -vzmin $MIN_VZ -vzmax $MAX_VZ
