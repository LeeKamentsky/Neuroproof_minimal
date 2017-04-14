# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (Algorithms)

set (SOURCES MergePriorityFunction.cpp)

add_library (Algorithms SHARED ${SOURCES})

# if (NOT BUILDEM_DIR STREQUAL "None")
#     add_dependencies (Classifier ${vigra_NAME} ${opencv_NAME})
# endif()
