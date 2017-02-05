SET(CDCVAM_EXAMPLES
   exampleCenterlineMesh
)



FOREACH(FILE ${CDCVAM_EXAMPLES})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} )
ENDFOREACH(FILE)





 





