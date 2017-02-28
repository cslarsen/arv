cmake_minimum_required(VERSION 3.5)
project(arv-cython C)

include(UseCython)

set_source_files_properties(${PROJECT_SOURCE_DIR}/arv.pyx
  PROPERTIES CYTHON_IS_CXX TRUE
)

cython_add_module(arv arv.pyx)

target_include_directories(arv
  PRIVATE ${PROJECT_SOURCE_DIR}/../cpp/
)

target_link_libraries(arv
  PUBLIC arv-cpp
)
