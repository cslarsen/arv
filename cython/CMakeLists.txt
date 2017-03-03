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

set_target_properties(_arv PROPERTIES
  VISIBILITY_INLINES_HIDDEN 1
  BUILD_WITH_INSTALL_RPATH FALSE
  INSTALL_RPATH "\$ORIGIN"
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS _arv LIBRARY DESTINATION lib/arv)
