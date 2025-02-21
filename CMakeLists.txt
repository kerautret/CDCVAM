PROJECT(CDCVAM)

cmake_minimum_required (VERSION 3.1)
cmake_policy(SET CMP0057 NEW) # the new interpretation of IN_LIST is by default not set (at least until  3.21.1).

SET(CDCVAM_SRC_COMP
 ${PROJECT_SOURCE_DIR}/src/GeodesicGraphComputer.cpp  
 ${PROJECT_SOURCE_DIR}/src/NormalAccumulator.cpp)



message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "Checking if DGtal is installed:")
FIND_PACKAGE(DGtal 0.9.3 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")



 # -----------------------------------------------------------------------------
 # CPP17
 # -----------------------------------------------------------------------------
set(CDCVAM_MIN_REQUIRED 17)
 if(NOT CMAKE_CXX_STANDARD)
   set(CMAKE_CXX_STANDARD ${CDCVAM_MIN_REQUIRED})
 else()
  # Throw if CMAKE_CXX_STANDARD is 98
  if(${CMAKE_CXX_STANDARD} EQUAL 98)
    message(FATAL_ERROR "CMAKE_CXX_STANDARD is set to ${CMAKE_CXX_STANDARD}, "
      "but DGtalTools requires at least ${CDCVAM_MIN_REQUIRED}.")
  endif()
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()






include_directories ("${PROJECT_SOURCE_DIR}/src")

# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )


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
