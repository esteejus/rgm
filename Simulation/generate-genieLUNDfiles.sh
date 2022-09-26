#!/bin/bash
#40Ar file names
#mgolden_G18_10a_02_11a_1000060120_598636.gst.root  mgolden_G18_10a_02_11a_1000180400_402962.gst.root  mgolden_G18_10a_02_11a_1000180400_598636.gst.root 
#C12 file names
#mgolden_G18_10a_02_11a_1000060120_207052.gst.root  mgolden_G18_10a_02_11a_1000060120_402962.gst.root  mgolden_G18_10a_02_11a_1000060120_598636.gst.root 

C12_PATH=/lustre19/expphy/volatile/clas12/mgolden/simulation/elements/12C/
Ar40_PATH=/lustre19/expphy/volatile/clas12/mgolden/simulation/elements/40Ar/

#10M in each gst file
#Usage of GENIE_to_LUND.C("inputFile","outputFile_prefix",#outputFiles,#eventsPerRun,"target",A,Z)

#root 'GENIE_to_LUND.C("${C12_PATH}/mgolden_G18_10a_02_11a_1000060120_207052.gst.root","c12_2gev",100,100000,"1-foil",6,6)'
#root 'GENIE_to_LUND.C("${C12_PATH}/mgolden_G18_10a_02_11a_1000060120_402962.gst.root","c12_4gev",100,100000,"1-foil",6,6)'
#root 'GENIE_to_LUND.C("${C12_PATH}/mgolden_G18_10a_02_11a_1000060120_598636.gst.root","c12_6gev",100,100000,"1-foil",6,6)'

#root 'GENIE_to_LUND.C("${Ar40_PATH}/40Ar/mgolden_G18_10a_02_11a_1000180400_402962.gst.root","ar40_4gev",100,100000,"Ar",20,20)'
#root 'GENIE_to_LUND.C("${Ar40_PATH}/mgolden_G18_10a_02_11a_1000180400_598636.gst.root","ar40_6gev",100,100000,"Ar",20,20)'
#root 'GENIE_to_LUND.C("${Ar40_PATH}/mgolden_G18_10a_02_11a_1000180400_207052.gst.root","ar40_2gev",100,100000,"Ar",20,20)'
