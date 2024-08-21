foreach(fnameSrc monitor_ep.cpp monitor_en.cpp monitor_epp.cpp monitor_epn.cpp monitorPID.cpp compare.cpp compare_epp.cpp compare_epn.cpp neff_h_epin_old.cpp neff_d_pcdn_old.cpp neff_d_pfdn_old.cpp skimmer.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks PkgConfig::hipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib NeutronVeto EventCut)
endforeach()

