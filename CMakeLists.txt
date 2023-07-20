cmake_minimum_required(VERSION 3.0)
project(RGM)
set (CLAS12ROOT $ENV{CLAS12ROOT})
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
set(LIBRARIES ${CMAKE_SOURCE_DIR}/libraries)


include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${INCLUDES} ${LIBRARIES})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

add_subdirectory(libraries)

add_subdirectory(Ana)
add_subdirectory(Analysis)
add_subdirectory(MagneticFieldStudy)
add_subdirectory(Monitoring)
add_subdirectory(NeutronVeto)
add_subdirectory(NeutronEfficiency)
