include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include})

## compile libraries in source directory
foreach(fname veto_functions.cpp)
	message(STATUS ${fname})
        string (REPLACE ".cpp" "" fnameLib ${fname})
        add_library(${fnameLib} ${fname})
        target_link_libraries(${fnameLib} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()
