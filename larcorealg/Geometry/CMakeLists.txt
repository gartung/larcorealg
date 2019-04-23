
cet_make( LIBRARIES
  larcorealg_CoreUtils
  
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  ${CETLIB}
  ${CETLIB_EXCEPT}
  ${CLHEP}
  ${ROOT_BASIC_LIB_LIST}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_GENVECTOR}
  )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
