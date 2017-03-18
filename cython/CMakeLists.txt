cmake_minimum_required(VERSION 3.5)
project(arv-cython-cpp CXX)

include(UseCython)

set_source_files_properties(${PROJECT_SOURCE_DIR}/_arv.pyx
  PROPERTIES CYTHON_IS_CXX TRUE
)

file(GLOB cpp_sources
  LIST_DIRECTORIES false
  ${CMAKE_SOURCE_DIR}/cpp/*.cpp
  ${CMAKE_SOURCE_DIR}/cpp/*.hpp
)

# Ready calling sparsehash configure
set(SHASH_BUILD_DIR "${CMAKE_BINARY_DIR}/cpp/sparsehash")
set(CONFIG_H        "${SHASH_BUILD_DIR}/src/config.h")
set(SPARSECONFIG_H  "${SHASH_BUILD_DIR}/internal/sparseconfig.h")

file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}")
file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}/src")
file(MAKE_DIRECTORY "${SHASH_BUILD_DIR}/internal")

# Sparsehash: Create config.h
add_custom_command(
  OUTPUT ${CONFIG_H}
  COMMAND ${CMAKE_SOURCE_DIR}/3rd-party/sparsehash/configure
  WORKING_DIRECTORY ${SHASH_BUILD_DIR}
  USES_TERMINAL
  COMMENT "Configuring Google Sparsehash"
)

# Sparsehash: Copy config.h to sparseconfig.h
add_custom_command(
  OUTPUT ${SPARSECONFIG_H}
  COMMAND ${CMAKE_COMMAND} -E copy ${CONFIG_H} ${SPARSECONFIG_H}
  MAIN_DEPENDENCY ${CONFIG_H}
)

cython_add_module(_arv _arv.pyx
  ${cpp_sources}
  ${SPARSECONFIG_H}
)

if (NOT APPLE)
  target_link_libraries(_arv PRIVATE ${PYTHON_LIBRARY})
endif()

target_include_directories(_arv
  PRIVATE ${CMAKE_SOURCE_DIR}/cpp # for cpp_sources
  PRIVATE ${CMAKE_BINARY_DIR}/cpp # sparseconfig.h
)

install(TARGETS _arv LIBRARY DESTINATION lib)
