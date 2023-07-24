foreach(fnameSrc many_plots.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameLib ${fnameSrc})
  add_library(${fnameLib} ${fnameSrc})
  target_link_libraries(${fnameLib} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()

foreach(fnameSrc Electron_Cuts.cpp Proton_Cuts.cpp ep_Kinematics.cpp epp_Kinematics.cpp  Andrew_skim.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib many_plots Clas12Ana)
endforeach()


