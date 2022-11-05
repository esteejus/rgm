cmake_minimum_required(VERSION 3.0)

project(Analysis)

set (CLAS12ROOT $ENV{CLAS12ROOT})
set(analysis_include ${CMAKE_SOURCE_DIR}/include )
set(analysis_libraries ${CMAKE_SOURCE_DIR}/libraries )


#include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${include})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include} ${analysis_libraries})

find_package(ROOT)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -O0 -ggdb")

add_subdirectory(source)

## compile libraries 
foreach(fnameSrc eventcut_old.cpp)
        message(STATUS ${fnameSRC})
        string (REPLACE ".cpp" "" fnameLib ${fnameSrc})
        add_library(${fnameLib} ${fnameSrc})
        target_link_libraries(${fnameLib} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()

## compile exacutables 
foreach(fnameSrc example_code.cpp tensorscalar.cpp) 
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib eventcut_old)
endforeach()


