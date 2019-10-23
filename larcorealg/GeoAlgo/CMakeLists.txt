cet_make(NO_DICTIONARY
         LIBRARIES
           ${ROOT_CORE}
           ${ROOT_PHYSICS}
         )

art_dictionary(DICTIONARY_LIBRARIES larcorealg_GeoAlgo)


install_headers()
install_source()
