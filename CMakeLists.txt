cmake_minimum_required(VERSION 3.0)
project(RGM)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set(include ${CMAKE_SOURCE_DIR}/Analysis/include )

include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${include})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

add_subdirectory(Analysis)

