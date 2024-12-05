foreach(fnameSrc Energy_Loss.cpp Energy_Loss_e.cpp Energy_Loss_CD.cpp Energy_Loss_FD.cpp Energy_Loss_FD_copyEdit.cpp Angular_Correction.cpp Angular_Correction_Hists.cpp Momentum_Correction.cpp Momentum_Correction_Hists.cpp Momentum_Correction_Hists_CD.cpp Momentum_Correction_Deuterium.cpp Efficiency_Resolution.cpp Get_Radiation_Simulation.cpp Efficiency_Resolution_Comp.cpp Resolution.cpp Resolution_Hist.cpp Resolution_Graph.cpp Vertex.cpp Vertex_Graph.cpp Q2.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} PkgConfig::hipo4 -lEG -lClas12Banks -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib Clas12Ana Clas12Debug reweighter many_plots)
endforeach()



