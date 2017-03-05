cmake_minimum_required(VERSION 3.5)
project(arv CXX)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Allow users to select which Python version to use
set(ARV_PYTHON_VERSION 3 CACHE
  STRING "Target Python version (e.g. 2, 3, 3.3, etc)")
set(Python_ADDITIONAL_VERSIONS ${ARV_PYTHON_VERSION})
find_package(PythonInterp)

# Python package
set(PACKAGE_VERSION "0.1")

# Set default build type to release
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Generate native code for release version
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND
    (CMAKE_COMPILER_IS_GNUCXX OR
     CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
     CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/cython")
add_subdirectory("${PROJECT_SOURCE_DIR}/python")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests")

message(STATUS "ARV_PYTHON_VERSION: ${ARV_PYTHON_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
message(STATUS "Python version: ${PYTHON_VERSION_STRING}")
