PROJECT(CDCVAM)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CDCVAM_SRC_COMP
 ${PROJECT_SOURCE_DIR}/src/GeodesicGraphComputer.cpp  
 ${PROJECT_SOURCE_DIR}/src/NormalAccumulator.cpp)

FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


include_directories ("${PROJECT_SOURCE_DIR}/src")

SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
      SET(CDCVAMLibDependencies ${CDCVAMLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )


if ( USE_PCL )
    find_package(PCL 1.3 REQUIRED COMPONENTS common segmentation)
    include_directories(${PCL_INCLUDE_DIRS})
    link_directories(${PCL_LIBRARY_DIRS})
    add_definitions(${PCL_DEFINITIONS})
    add_definitions(-DUSE_PCL=${USE_PCL})
    SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES})
else ()
    SET(PCLLib "")
endif()


add_subdirectory(bin)
add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(tools)
