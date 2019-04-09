# define ROOTLIB

set( ROOTLIB ROOT::Core ROOT::RIO ROOT::Net ROOT::Hist ROOT::Graf ROOT::Graf3d ROOT::Gpad ROOT::Tree ROOT::Rint ROOT::Postscript ROOT::Matrix ROOT::Physics ROOT::MathCore ROOT::Thread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} ${ROOTLIB} ${ROOT_LIBRARIES})

cet_make( LIBRARIES NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES larcorealg_GeoAlgo )


install_headers()
install_source()

