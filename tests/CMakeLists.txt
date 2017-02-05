SET(CDCVAM_TESTS
   testCenterLineHelper
   testGeodesicGraph
   testGeodesicGraphVertexRepr
   testNormalAccumulator
   testNormalAccumulatorNormals
)

if (  WITH_VISU3D_QGLVIEWER )
   SET(CDCVAM_TESTS ${CDCVAM_TESTS} testGeodesicGraphDisplay)
endif ( WITH_VISU3D_QGLVIEWER )  



if ( USE_PCL )
   SET( CDCVAM_TESTS_PCM
    ${CDCVAM_TESTS_PCL}
    # ...
)
endif()


FOREACH(FILE ${CDCVAM_TESTS})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} )
ENDFOREACH(FILE)

FOREACH(FILE ${CDCVAM_TESTS_PCL})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





