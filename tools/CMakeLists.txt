SET(CDCVAM_TOOLS
 statisticsCC
 statisticsCentricity  
 statisticsRadiusErrDistr 
 statisticsRadiusErrEvol 
 statisticsRadiusEstim
)



FOREACH(FILE ${CDCVAM_TOOLS})
  add_executable(${FILE} ${FILE} ${CDCVAM_SRC_COMP})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${CDCVAMLibDependencies} ${PCLLib})
ENDFOREACH(FILE)





 





