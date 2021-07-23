#!/bin/bash
SOURCE=$1
DIR=$2
RUN1=$3
RUN2=$4

for RUN in $(seq $RUN1 $RUN2)
do
    echo "${SOURCE} ${RUN}"
    clas12root "CLAS12Skimmer.C(\"${SOURCE}${RUN}*\",\"${DIR}${RUN}_skim.hipo\")"
done

