cet_make(
  SUBDIRS
    "details"
  LIBRARIES
          ${MF_MESSAGELOGGER}
          ${FHICLCPP}
          cetlib
          cetlib_except
          ROOT::Core
          ROOT::Physics
          ROOT::Geom
          ROOT::GenVector
          )

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")

add_subdirectory(gdml)
