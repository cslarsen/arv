cmake_minimum_required(VERSION 3.5)
project(arv-cpp CXX)

file(GLOB sources
  LIST_DIRECTORIES false
  ${PROJECT_SOURCE_DIR}/*.cpp
  ${PROJECT_SOURCE_DIR}/*.hpp
)

set(SHASH_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/sparsehash")
set(CONFIG_H        "${SHASH_BUILD_DIR}/src/config.h")
set(SPARSECONFIG_H  "${SHASH_BUILD_DIR}/internal/sparseconfig.h")

file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}")
file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}/src")
file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}/internal")

# Create config.h
add_custom_command(
  OUTPUT ${CONFIG_H}
  COMMAND ${CMAKE_SOURCE_DIR}/3rd-party/sparsehash/configure
  WORKING_DIRECTORY ${SHASH_BUILD_DIR}
  USES_TERMINAL
  COMMENT "Configuring Google Sparsehash"
)

# Copy config.h to sparseconfig.h
add_custom_command(
  OUTPUT ${SPARSECONFIG_H}
  COMMAND ${CMAKE_COMMAND} -E copy ${CONFIG_H} ${SPARSECONFIG_H}
  MAIN_DEPENDENCY ${CONFIG_H}
)

add_library(arv SHARED
  ${sources}
  ${SPARSECONFIG_H}
)

target_include_directories(arv
  PUBLIC ${PROJECT_SOURCE_DIR}
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)

# Hide all non-public symbols
set_target_properties(arv PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1
)

install(TARGETS arv LIBRARY DESTINATION lib/arv)
