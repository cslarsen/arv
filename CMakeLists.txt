cmake_minimum_required(VERSION 3.5)
project(arv CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

add_subdirectory("${PROJECT_SOURCE_DIR}/cpp")
add_subdirectory("${PROJECT_SOURCE_DIR}/cython")
