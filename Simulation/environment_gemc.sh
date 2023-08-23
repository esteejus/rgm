#!/bin/bash                                                                                    
source /group/clas12/packages/setup.sh

#set our Environment on JLab Farm
module purge
module load cmake
module load sqlite/dev
module load clas12/pro
module switch coatjava/10.0.1
module switch gemc/5.2
module switch root/6.20.04 #needed to work with GCF compilaiton 
