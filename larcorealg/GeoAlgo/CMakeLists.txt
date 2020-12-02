cet_make(NO_DICTIONARY ### MIGRATE-ACTION-RECOMMENDED: use cet_make_library(), build_dictionary(), basic_plugin() with explicit source lists
         LIBRARIES
           ROOT::Core
           ROOT::Physics
         )

art_dictionary(DICTIONARY_LIBRARIES larcorealg_GeoAlgo)


install_headers()
install_source()
