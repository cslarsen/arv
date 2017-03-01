cmake_minimum_required(VERSION 3.5)
project(arv-cpp CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(google_densehash REQUIRED)

file(GLOB sources
  LIST_DIRECTORIES false
  ${PROJECT_SOURCE_DIR}/*.cpp
  ${PROJECT_SOURCE_DIR}/*.hpp
)

add_library(arv-cpp SHARED ${sources})

target_include_directories(arv-cpp
  PUBLIC ${PROJECT_SOURCE_DIR}
  PRIVATE ${google_densehash_INCLUDE_DIR}
)

set_target_properties(arv-cpp PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1
)
