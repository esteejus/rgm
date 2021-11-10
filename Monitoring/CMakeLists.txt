cmake_minimum_required(VERSION 3.0)
project(monitorPID)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set(analysis_include ${CMAKE_SOURCE_DIR}/include )

include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

add_subdirectory(source)

# compile executables in source directory
foreach(fnameSrc monitor_ep.cpp monitor_epp.cpp monitor_epn.cpp compare_epp.cpp compare_epn.cpp skimmer.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib eventcut)
endforeach()

#Other exacutables that haven't been pushed
#foreach(fnameSrc monitorPID.cpp compare_gen.cpp compare.cpp monitor_6GeV_H.cpp monitor_6GeV_D.cpp monitor_noFTOF.cpp pmiss.cpp)
