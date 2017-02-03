SET(CDCVAM_BIN
 centerLineGeodesicGraph
 compAccFromMesh
 displayAccuMesh
)

if (  WITH_VISU3D_QGLVIEWER )
   SET(CDCVAM_BIN ${CDCVAM_BIN} illustrationGraphAllSteps)
endif ( WITH_VISU3D_QGLVIEWER )  


if ( USE_PCL )
   SET( CDCVAM_BIN
    ${CDCVAM_BIN}
    compAccFromSDP )
endif()


FOREACH(FILE ${CDCVAM_BIN})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





