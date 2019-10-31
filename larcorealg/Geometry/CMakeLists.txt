cet_make(
  SUBDIRS
    "details"
  LIBRARIES
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

# NOTE differently from the usual pattern, `larcorealg_Geometry_PixelPlane`
# library depends on `larcorealg_Geometry`
# (the non-disruptive alternative being move it up one directory, with the
# drawback that core geometry material get spread around; or to dump everything
# in `Geometry`, which is large enough already)
add_subdirectory(PixelPlane)

add_subdirectory(gdml)
