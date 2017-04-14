# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (DataStructures)

set (SOURCES Stack.cpp StackLearn.cpp StackPredict)

add_library (DataStructures SHARED ${SOURCES})

target_link_libraries (DataStructures Algorithms)

# if (NOT BUILDEM_DIR STREQUAL "None")
#     add_dependencies (Classifier ${vigra_NAME} ${opencv_NAME})
# endif()
