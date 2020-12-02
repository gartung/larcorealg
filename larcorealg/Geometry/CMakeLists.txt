cet_make_library(
  LIBRARY_NAME larcorealg_Geometry
  SOURCE
    AuxDetChannelMapAlg.cxx
    AuxDetGeo.cxx
    AuxDetGeometryCore.cxx
    AuxDetSensitiveGeo.cxx
    BoxBoundedGeo.cxx
    ChannelMapAlg.cxx
    ChannelMapStandardAlg.cxx
    CryostatGeo.cxx
    DriftPartitions.cxx
    GeometryBuilderStandard.cxx
    GeometryBuilderWireless.cxx
    GeometryCore.cxx
    GeoNodePath.cxx
    GeoObjectSorter.cxx
    GeoObjectSorterStandard.cxx
    GeoVectorLocalTransformation.cxx
    LocalTransformation.cxx
    OpDetGeo.cxx
    PlaneGeo.cxx
    StandaloneGeometrySetup.cxx
    TPCGeo.cxx
    WireGeo.cxx
  LIBRARIES
          larcoreobj::larcoreobj_SimpleTypesAndConstants
          messagefacility::MF_MessageLogger
          fhiclcpp::fhiclcpp
          cetlib::cetlib
          cetlib_except::cetlib_except
          ROOT::Core
          ROOT::Physics
          ROOT::Geom
          ROOT::GenVector
          CLHEP::CLHEP
        )

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
install_source(SUBDIRS "details")

add_subdirectory(gdml)
