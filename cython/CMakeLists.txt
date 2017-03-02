cmake_minimum_required(VERSION 3.5)
project(arv-cython CXX)

include(UseCython)

set_source_files_properties(${PROJECT_SOURCE_DIR}/_arv.pyx
  PROPERTIES CYTHON_IS_CXX TRUE
)

cython_add_module(_arv _arv.pyx)

set(libs PRIVATE arv)

if (NOT APPLE)
  list(APPEND libs PRIVATE ${PYTHON_LIBRARY})
endif()

target_link_libraries(_arv ${libs})

install(TARGETS _arv LIBRARY DESTINATION lib/arv)
