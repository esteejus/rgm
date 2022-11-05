
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${CLAS12ROOT}/hipo4 ${CLAS12ROOT}/rcdb/cpp/include ${analysis_include})

## compile libraries in source directory
foreach(fnameSrc eventcut.cpp)
	message(STATUS ${fnameSRC})
        string (REPLACE ".cpp" "" fnameLib ${fnameSrc})
        add_library(${fnameLib} ${fnameSrc})
        target_link_libraries(${fnameLib} ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
endforeach()

# Add library for target info
#add_library(eventcut eventcut.cpp)
#target_link_libraries(eventcut ${ROOT_LIBRARIES} -lEG -lClas12Banks -lHipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib)
