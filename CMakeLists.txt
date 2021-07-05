PROJECT(CDCVAM)

cmake_minimum_required(VERSION 3.19)

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
 # CPP11
 # -----------------------------------------------------------------------------
 set (CMAKE_CXX_STANDARD 11)
 set (CMAKE_CXX_STANDARD_REQUIRED TRUE)





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
