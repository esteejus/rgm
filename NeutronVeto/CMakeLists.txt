cmake_minimum_required(VERSION 3.0)

project(neutronveto)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set(analysis_include ${CMAKE_SOURCE_DIR}/include )
set(analysis_libraries ${CMAKE_SOURCE_DIR}/libraries )

include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include} ${analysis_libraries})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

#add_subdirectory(source)

foreach (fname D_getfeatures.cpp)
  message(STATUS ${fname})
  string (REPLACE ".cpp" "" fnameExe ${fname})
  add_executable(${fnameExe} ${fname})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()
