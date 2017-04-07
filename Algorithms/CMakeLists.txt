# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (Algorithms)

# IF(NOT CMAKE_BUILD_TYPE)
#     SET (CMAKE_BUILD_TYPE Release)
# ENDIF(NOT CMAKE_BUILD_TYPE)
# 
# SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
# SET(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
# SET(CMAKE_CXX_LINK_FLAGS "-O3")
# SET(CMAKE_DEBUG_POSTFIX "-g")
# 

   #set (CMAKE_C_COMPILER ${CMAKE_SOURCE_DIR}/cilkplus-install/bin/gcc) 
   #set (CMAKE_CXX_COMPILER ${CMAKE_SOURCE_DIR}/cilkplus-install/bin/g++)
 #SET(CMAKE_CXX_FLAGS "-fcilkplus")
 #SET(CMAKE_CXX_LINK_FLAGS "-lcilkrts")
   #SET(CMAKE_CXX_FLAGS "-fcilkplus -O3")
   #SET(CMAKE_CXX_LINK_FLAGS "-lcilkrts")
   #set (CMAKE_C_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/../cilkplus-install/bin/gcc") 
   #set (CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/../cilkplus-install/bin/g++")
   SET(CMAKE_CXX_FLAGS "-fcilkplus -O3")
   SET(CMAKE_CXX_LINK_FLAGS "-lcilkrts")
set (SOURCES MergePriorityFunction.cpp)

add_library (Algorithms SHARED ${SOURCES})

# if (NOT BUILDEM_DIR STREQUAL "None")
#     add_dependencies (Classifier ${vigra_NAME} ${opencv_NAME})
# endif()
