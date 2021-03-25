#!/bin/bash

cd /work/clas12/users/esteejus/clas12analysis/ForFarm/LowEnergyDeuteron

gemc -USE_GUI=0 -N=10 -INPUT_GEN_FILE="lund, ./lundfiles/lund_qe_1.dat" -OUTPUT="evio, ./eviofiles/out_1.ev" /group/clas12/gemc/4.4.1/config/rgb_spring2019.gcard
evio2hipo -t -1 -s -1 -i ./eviofiles/out_1.ev -o ./mchipo/qe_1.hipo
recon-util -y /group/clas12/gemc/4.4.1/config/rgb_spring2019.yaml -i ./mchipo/qe_1.hipo -o ./reconhipo/recon_qe_1.hipo





#gemc -USE_GUI=0 -N=10000 -INPUT_GEN_FILE="lund, ./lundfiles/lund_qe_1.dat" /group/clas12/gemc/4.4.1/config/rgb_spring2019.gcard
#evio2hipo -t -1 -s -1 -i out.ev -o qe_1.hipo
#recon-util -y /group/clas12/gemc/4.4.1/config/rgb_spring2019.yaml -i qe_1.hipo -o recon_qe_1.hipo
