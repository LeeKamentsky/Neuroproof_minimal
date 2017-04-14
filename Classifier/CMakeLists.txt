# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (Classifier)

set (SOURCES opencvRFclassifier.cpp vigraRFclassifier.cpp)

add_library (Classifier SHARED ${SOURCES})

# if (NOT BUILDEM_DIR STREQUAL "None")
#     add_dependencies (Classifier ${vigra_NAME} ${opencv_NAME})
# endif()
