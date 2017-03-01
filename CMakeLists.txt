cmake_minimum_required(VERSION 3.5)
project(arv CXX)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Allow users to select which Python version to use
set(ARV_PYTHON_VERSION 3 CACHE STRING "Target Python version (e.g. 2 or 3)")
message(STATUS "ARV_PYTHON_VERSION ${ARV_PYTHON_VERSION}")
set(Python_ADDITIONAL_VERSIONS ${ARV_PYTHON_VERSION})

# Generate native code for release version
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND
    (CMAKE_COMPILER_IS_GNUCXX OR
     CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
     CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
  message(STATUS "Adding -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/cpp")
add_subdirectory("${PROJECT_SOURCE_DIR}/cython")
