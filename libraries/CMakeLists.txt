# add library names
add_library(NeutronVeto neutron-veto/veto_functions.cpp)
add_library(EventCut eventcut/eventcut.cpp)
target_link_libraries(NeutronVeto ${ROOT_LIBRARIES})
target_link_libraries(EventCut ${ROOT_LIBRARIES})

add_library(Efficiency efficiency/efficiency.cpp)
target_link_libraries(Efficiency ${ROOT_LIBRARIES})

add_library(Clas12Ana clas12ana/clas12ana.cpp)
add_library(Clas12Debug clas12debug/clas12debug.cpp)
target_link_libraries(Clas12Ana ${ROOT_LIBRARIES})
target_link_libraries(Clas12Debug ${ROOT_LIBRARIES})

add_library(eNCrossSection simulation_reweighting/eNCrossSection.cc)
target_link_libraries(eNCrossSection ${ROOT_LIBRARIES})
add_library(gcfSRC simulation_reweighting/gcfSRC.cc)
target_link_libraries(gcfSRC ${ROOT_LIBRARIES})
add_library(reweighter simulation_reweighting/reweighter.cpp)
target_link_libraries(reweighter gcfSRC eNCrossSection ${ROOT_LIBRARIES}  -lEG -lClas12Banks PkgConfig::hipo4 -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib Clas12Ana Clas12Debug NeutronVeto)

