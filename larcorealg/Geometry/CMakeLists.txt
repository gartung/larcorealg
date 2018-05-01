
cet_make( LIBRARIES
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except
  ${CLHEP}
  ${ROOT_BASIC_LIB_LIST}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_EG}
  ${ROOT_GENVECTOR}
  )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
