cmake_minimum_required(VERSION 3.0)
project(monitorPID)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set (exename monitorPID)

include_directories(${CLAS12ROOT}/Clas12Banks)
include_directories(${CLAS12ROOT}/Clas12Root)
include_directories(${CLAS12ROOT}/hipo4)
include_directories(${CLAS12ROOT}/rcdb/cpp/include)


find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


add_executable(${exename}  monitorPID.cpp )
 
target_link_libraries(${exename} ${ROOT_LIBRARIES} -lEG)
target_link_libraries(${exename} -lClas12Banks -lHipo4 -lClas12Root)
target_link_libraries(${exename} -L${CLAS12ROOT}/lib)
target_link_libraries(${exename} -L${CLAS12ROOT}/ccdb/lib)