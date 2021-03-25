cmake_minimum_required(VERSION 3.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name and version
project(Monitoring VERSION 1.0)

find_package(ROOT REQUIRED)
find_package(CLAS12ROOT REQUIRED HINTS /w/hallb-scifs17exp/clas12/users/esteejus/clas12root/)

if(DEFINED ENV{CLAS12ROOT}) 
 set(CLAS12ROOT_DIR $ENV{CLAS12ROOT})
endif()

include_directories(${CLAS12ROOT_DIR}/Clas12Root)
include_directories(${CLAS12ROOT_DIR}/hipo4)
include_directories(${CLAS12ROOT_DIR}/Clas12Banks)
include_directories(${CLAS12ROOT_DIR}/rcdb/cpp/include/)

link_directories(${CLAS12ROOT_DIR}/Clas12Root)
link_directories(${CLAS12ROOT_DIR}/Clas12Banks)

add_library(hipochain ${CLAS12ROOT_DIR}/Clas12Root/HipoChain.h ${CLAS12ROOT_DIR}/Clas12Root/HipoChain.cpp)

add_library(clas12reader ${CLAS12ROOT_DIR}/Clas12Banks/clas12reader.h  ${CLAS12ROOT_DIR}/Clas12Banks/clas12reader.cpp) 

target_link_libraries(clas12reader 
ROOT::Core
ROOT::Gpad
ROOT::Graf3d
ROOT::Graf
ROOT::Hist
ROOT::Imt
ROOT::MathCore
ROOT::Matrix
ROOT::MultiProc
ROOT::Net
ROOT::Physics
ROOT::Postscript
ROOT::RIO
ROOT::ROOTDataFrame
ROOT::ROOTVecOps
ROOT::Rint
ROOT::EG
ROOT::Thread
ROOT::TreePlayer
ROOT::Tree)

target_link_libraries(hipochain 
ROOT::Core
ROOT::Gpad
ROOT::Graf3d
ROOT::Graf
ROOT::Hist
ROOT::Imt
ROOT::MathCore
ROOT::Matrix
ROOT::MultiProc
ROOT::Net
ROOT::Physics
ROOT::Postscript
ROOT::EG
ROOT::RIO
ROOT::ROOTDataFrame
ROOT::ROOTVecOps
ROOT::Rint
ROOT::Thread
ROOT::TreePlayer
ROOT::Tree
)

# add the executable
add_executable(Mon LowEnergyReader.cpp)
target_link_libraries(Mon clas12reader hipochain 
ROOT::Core
ROOT::Gpad
ROOT::Graf3d
ROOT::Graf
ROOT::Hist
ROOT::Imt
ROOT::MathCore
ROOT::Matrix
ROOT::MultiProc
ROOT::Net
ROOT::Physics
ROOT::Postscript
ROOT::EG
ROOT::RIO
ROOT::ROOTDataFrame
ROOT::ROOTVecOps
ROOT::Rint
ROOT::Thread
ROOT::TreePlayer
ROOT::Tree
)
