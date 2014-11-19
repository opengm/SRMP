# default build-type (release)
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

# check cmake version
cmake_minimum_required(VERSION 2.8)

# project name
project (SRMP)

# set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# set source files
SET(SRMP_SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/AddEdges.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/CMP.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/debug.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/energy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/MPLP.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/order.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/SRMP.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Algs/triplets.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/GeneralType.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/PairwiseDualType.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/PairwiseType.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/PatternType.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/PottsType.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FactorTypes/SharedPairwiseType.cpp
)

# create library
add_library (srmp ${SRMP_SRC_FILES})

# point to include directory
target_include_directories (srmp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# build example
add_executable (srmpExample ${CMAKE_CURRENT_SOURCE_DIR}/src/examples/example.cpp)
# point to include directory
target_include_directories (srmpExample PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# link srmp library
target_link_libraries (srmpExample LINK_PUBLIC srmp)

#link rt library
if(UNIX AND NOT APPLE)
   find_library(RT rt)
   target_link_libraries (srmpExample LINK_PUBLIC rt)
endif(UNIX AND NOT APPLE)