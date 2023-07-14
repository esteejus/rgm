#!/bin/csh                                                                                    
source /group/clas12/packages/setup.csh
#set our Environment on JLab Farm
module purge
module load cmake
module load sqlite/dev
module load clas12/pro
module switch coatjava/10.0.2
