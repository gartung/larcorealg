
cet_make( LIBRARIES
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except
  ${CLHEP_LIBRARIES}
  ${ROOT_LIBRARIES}
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ROOT::EG
  ROOT::GenVector
  )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
