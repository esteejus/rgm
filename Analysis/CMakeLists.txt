foreach(fnameSrc example_code.cpp tensorscalar.cpp tensorscalar_mc.cpp) 
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks PkgConfig::hipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib -lTMVA NeutronVeto EventCut)
endforeach()
