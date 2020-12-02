cet_make_library(
         SOURCE
           GeoAABox.cxx
           GeoAlgo.cxx
           GeoCone.cxx
           GeoCylinder.cxx
           GeoDirectedLine.cxx
           GeoHalfLine.cxx
           GeoLine.cxx
           GeoLineSegment.cxx
           GeoObjCollection.cxx
           GeoSphere.cxx
           GeoTrajectory.cxx
           GeoVector.cxx
         LIBRARIES
           ROOT::Core
           ROOT::Physics
         )

art_dictionary(DICTIONARY_LIBRARIES larcorealg_GeoAlgo)


install_headers()
install_source()
