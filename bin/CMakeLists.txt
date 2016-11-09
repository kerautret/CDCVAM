SET(CDCVAM_BIN
 centerLineGeodesicGraph
 compAccFromMesh
 illustrationGraphAllSteps
 displayAccuMesh
)



if ( USE_PCL )
   SET( CDCVAM_BIN
    ${CDCVAM_BIN}
    compAccFromSDP )
endif()


FOREACH(FILE ${CDCVAM_BIN})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





