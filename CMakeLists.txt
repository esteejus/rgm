cmake_minimum_required(VERSION 3.29)
project(RGM)
set(CMAKE_CXX_STANDARD 17)

set(CLAS12ROOT $ENV{CLAS12ROOT})
set(HIPO $ENV{HIPO})
find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")
find_package(PkgConfig REQUIRED)
pkg_check_modules(hipo4 REQUIRED IMPORTED_TARGET hipo4)


set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/libraries/simulation_reweighting ${HIPO}/include)
set(LIBRARIES ${CMAKE_SOURCE_DIR}/libraries ${HIPO}/lib)
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})


add_definitions(-D_CLAS12ANA_DIR="${CMAKE_SOURCE_DIR}")

add_subdirectory(libraries)
add_subdirectory(Ana)
add_subdirectory(Analysis)
add_subdirectory(Monitoring)
add_subdirectory(NeutronVeto)
#add_subdirectory(NeutronEfficiency)
