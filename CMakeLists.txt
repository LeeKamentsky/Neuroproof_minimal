CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
   set (CMAKE_C_COMPILER ${CMAKE_SOURCE_DIR}/cilkplus-install/bin/gcc) 
   set (CMAKE_CXX_COMPILER ${CMAKE_SOURCE_DIR}/cilkplus-install/bin/g++)
   link_directories (${CMAKE_SOURCE_DIR}/cilkplus-install/lib64)
project (Neuroproof_minimal)

#  ask to use -DCMAKE_PREFIX_PATH=/groups/scheffer/home/paragt/miniconda/envs/tst_my/

#   include_directories (/groups/scheffer/home/paragt/miniconda/envs/tst_my/include/)

   include_directories (${CMAKE_PREFIX_PATH}/include)
   #include_directories (${CMAKE_PREFIX_PATH}/include)
   include_directories (Algorithms)
   include_directories (Utilities)
   include_directories (Watershed)
   include_directories (Classifier)
   include_directories (DataStructures)
   include_directories (FeatureManager)

#    link_directories (/groups/scheffer/home/paragt/miniconda/envs/tst_my/lib) 
   link_directories (${CMAKE_SOURCE_DIR}/cilkplus-install/lib64)
   link_directories (${CMAKE_PREFIX_PATH}/lib)
   #link_directories (/home/armafire/tools/vigra-install/lib64)
   #link_directories (/usr/lib/x86_64-linux-gnu)
   #link_directories (/lib/x86_64-linux-gnu)
   add_subdirectory (Algorithms)
   add_subdirectory (Classifier)
   add_subdirectory (DataStructures)
   add_subdirectory (FeatureManager)
 
   #set (CMAKE_C_COMPILER "/home/armafire/tools/cilkplus-install/bin/cc")
   #set (CMAKE_CXX_COMPILER "/home/armafire/tools/cilkplus-install/bin/c++")
   #set(CMAKE_CXX_FLAGS "-fcilkplus -std=c++0x")
   #set(CMAKE_CXX_FLAGS "-fcilkplus -lcilkrts -std=gnu++98 -g -mfpmath=sse -msse2")
   set(CMAKE_CXX_FLAGS "-fcilkplus -std=gnu++98 -O3")
   set(CMAKE_CXX_LINK_FLAGS "-lcilkrts")
   #set(CMAKE_CXX_FLAGS "-fcilkplus -std=c++11")
#    set(CMAKE_CXX_STANDARD 11)
#    set(CMAKE_CXX_STANDARD_REQUIRED ON)


#    add_executable (predict_pixels predict_pixels.cpp pixel_detector.cpp) 
#    add_executable (test_vigra_watershed test_vigra_watershed.cpp) 


#    add_library(Classifier SHARED Classifier/opencvRFclassifier.cpp
   add_executable (NeuroProof_stack NeuroProof_stack.cpp)  
   add_executable (NeuroProof_stack_learn NeuroProof_stack_learn.cpp)  

#libvigraimpex.so.5
   target_link_libraries (NeuroProof_stack Algorithms Classifier DataStructures FeatureManager vigraimpex opencv_ml opencv_core boost_python boost_thread boost_system boost_chrono python2.7 png hdf5 hdf5_hl cilkrts)  
   target_link_libraries (NeuroProof_stack_learn Algorithms Classifier DataStructures FeatureManager vigraimpex opencv_ml opencv_core boost_python boost_thread boost_system boost_chrono python2.7 png hdf5 hdf5_hl cilkrts)  
#    target_link_libraries (predict_pixels Classifier SemiSupervised Utilities vigraimpex opencv_ml opencv_core boost_thread boost_system boost_chrono hdf5 hdf5_hl)  
#    target_link_libraries (test_vigra_watershed Utilities vigraimpex opencv_ml opencv_core boost_thread boost_system boost_chrono hdf5 hdf5_hl)  


############################################################################
# endif()
