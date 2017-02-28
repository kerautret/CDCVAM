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
   SET( CDCVAM_TESTS
    ${CDCVAM_TESTS}
    # ...
)
endif()




FOREACH(FILE ${CDCVAM_TESTS})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  if ( USE_PCL )           
    target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} ${PCLLib})
  else ()
    target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} )
  endif( USE_PCL )  
ENDFOREACH(FILE)








 





