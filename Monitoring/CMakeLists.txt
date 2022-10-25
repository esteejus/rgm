cmake_minimum_required(VERSION 3.0)
project(monitorPID)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set(analysis_include ${CMAKE_SOURCE_DIR}/include )
set(analysis_libraries ${CMAKE_SOURCE_DIR}/libraries )

include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include} ${analysis_libraries})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

add_subdirectory(source)
#add_subdirectory(Andrew_Analysis)
#add_subdirectory(Neutrons_CD)
#add_subdirectory(Neutrons_ECAL)
#add_subdirectory(Protons)
#add_subdirectory(Q2_dep)
#add_subdirectory(GCF_tools)

#compile executables in source directory
foreach(fnameSrc monitor_ep.cpp monitor_epp.cpp monitor_epn.cpp monitorPID.cpp compare.cpp compare_epp.cpp compare_epn.cpp neff_h_epin.cpp neff_d_pcdn.cpp neff_d_pfdn.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib eventcut)
endforeach()

