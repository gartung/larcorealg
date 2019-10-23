cet_make(LIBRARIES
          larcorealg_CoreUtils
          ${MF_MESSAGELOGGER}
          ${FHICLCPP}
          ${CETLIB}
          ${CETLIB_EXCEPT}
          ${ROOT_CORE}
          ${ROOT_PHYSICS}
          ${ROOT_GEOM}
          ${ROOT_GENVECTOR}
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
