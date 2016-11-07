SET(CDCVAM_TESTS
   testCenterLineHelper
   testGeodesicGraph
   testGeodesicGraphDisplay
   testGeodesicGraphVertexRepr
   testNormalAccumulator
   testNormalAccumulatorNormals
)



if ( USE_PCL )
   SET( CDCVAM_TESTS
    ${CDCVAM_TESTS}
    # ...
)
endif()


FOREACH(FILE ${CDCVAM_TESTS})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





