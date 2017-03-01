cmake_minimum_required(VERSION 3.5)
project(arv CXX)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Adding -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/cpp")
add_subdirectory("${PROJECT_SOURCE_DIR}/cython")
