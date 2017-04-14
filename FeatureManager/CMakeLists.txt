# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (FeatureManaager)

set (SOURCES Features.cpp FeatureManager.cpp)

add_library (FeatureManager SHARED ${SOURCES})
