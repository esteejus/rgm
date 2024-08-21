project(ProtonEfficiency)
# compile executables in source directory
foreach(fnameSrc ProtonEfficiency.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks PkgConfig::hipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib Clas12Ana Clas12Debug -lTMVA NeutronVeto EventCut)
endforeach()

