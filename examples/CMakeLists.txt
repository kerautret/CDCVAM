SET(CDCVAM_EXAMPLES
   exampleCenterlineMesh
)





FOREACH(FILE ${CDCVAM_EXAMPLES})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  if ( USE_PCL )
    target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} ${PCLLib})
  else ()
      target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} )
  endif( USE_PCL)
ENDFOREACH(FILE)





 





