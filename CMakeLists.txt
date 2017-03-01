cmake_minimum_required(VERSION 3.5)
project(arv CXX)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Allow users to select which Python version to use
set(ARV_PYTHON_VERSION 3 CACHE STRING "Target Python version (e.g. 2 or 3)")
set(Python_ADDITIONAL_VERSIONS ${ARV_PYTHON_VERSION})

# Generate native code for release version
if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Adding -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/cpp")
add_subdirectory("${PROJECT_SOURCE_DIR}/cython")
