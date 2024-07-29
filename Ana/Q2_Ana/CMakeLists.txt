add_library(many_plots many_plots.cpp)
target_link_libraries(many_plots ${ROOT_LIBRARIES})

foreach(fnameSrc Fast_Plots.cpp Andrew_skim.cpp test_skim.cpp Central_Detector_Proton.cpp Q2_dependence.cpp rew_guas.cpp Deuterium_Lowmmiss_Inspection.cpp Proton_Resolution.cpp MC_Debug.cpp p_LUND.cpp Electron_DeltaOmega.cpp skim_elastic.cpp Electron_Resolution.cpp  Electron_DeltaE.cpp Deuterium_MissingMass_Bins.cpp PID_check.cpp skim_deepFDpCD.cpp Momentum_Corrections.cpp Energy_Loss.cpp iso_p_LUND.cpp Energy_Loss_FD.cpp Energy_Loss_CD.cpp Energy_Loss_e.cpp Momentum_Corrections_epi.cpp Momentum_Corrections_epi_hist.cpp skim_epi.cpp Analysis_Note_Q2.cpp Proton_Fid.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} PkgConfig::hipo4 -lEG -lClas12Banks -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib Clas12Ana Clas12Debug reweighter many_plots)
endforeach()
