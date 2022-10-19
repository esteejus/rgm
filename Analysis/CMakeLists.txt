
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${include})

## compile libraries 
foreach(fnameSrc eventcut_old.cpp)
        message(STATUS ${fnameSRC})
        string (REPLACE ".cpp" "" fnameLib ${fnameSrc})
        add_library(${fnameLib} ${fnameSrc})
        target_link_libraries(${fnameLib} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()

## compile exacutables 
foreach(fnameSrc example_code.cpp) 
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib eventcut_old)
endforeach()


