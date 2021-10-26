#!/bin/csh                                                                                    
source /group/clas12/packages/setup.csh
#set our Environment on JLab Farm
module purge
module load cmake
module load sqlite/4.4.1
module load clas12/pro
module switch clas12root/1.7.4
