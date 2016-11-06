SET(TUBE_ANALYSE_BIN
 centerLineFMMc
 centerLineGeodesicGraph
 compAccFromMesh
# compAccFromVol
 compBallFromAccMesh 
 illustrationGraphAllSteps
 displayAccuMesh
 statisticsCC
 statisticsCentricity  
 statisticsRadiusErrDistr 
 statisticsRadiusErrEvol 
 statisticsRadiusEstim
)



if ( USE_PCL )
   SET( TUBE_ANALYSE_BIN
    ${TUBE_ANALYSE_BIN}
    compAccFromSDP )
endif()


FOREACH(FILE ${TUBE_ANALYSE_BIN})
  add_executable(${FILE} ${FILE} ${TUBE_ANALYSE_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





