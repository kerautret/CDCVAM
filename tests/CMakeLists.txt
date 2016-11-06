SET(TUBE_ANALYSE_TESTS
   testCenterLineHelper
   testGeodesicGraph
   testGeodesicGraphDisplay
   testGeodesicGraphVertexRepr
   testNormalAccumulator
   testNormalAccumulatorNormals
)



if ( USE_PCL )
   SET( TUBE_ANALYSE_TESTS
    ${TUBE_ANALYSE_TESTS}
    # ...
)
endif()


FOREACH(FILE ${TUBE_ANALYSE_TESTS})
  add_executable(${FILE} ${FILE} ${TUBE_ANALYSE_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





