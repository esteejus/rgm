# add library names
add_library(NeutronVeto neutron-veto/veto_functions.cpp)
add_library(EventCut eventcut/eventcut.cpp)
target_link_libraries(NeutronVeto ${ROOT_LIBRARIES})
target_link_libraries(EventCut ${ROOT_LIBRARIES})

add_library(Clas12Ana clas12ana/clas12ana.cpp)
target_link_libraries(Clas12Ana ${ROOT_LIBRARIES})

