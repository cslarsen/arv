cmake_minimum_required(VERSION 3.5)
project(arv-python CXX)

configure_file(arv/__init__.py arv/__init__.py COPYONLY) # only for testing
configure_file(setup.py.in setup.py)

file(GLOB python_files
  ${CMAKE_CURRENT_SOURCE_DIR}/arv/*.py
  ${CMAKE_CURRENT_BUILD_DIR}/setup.py
)

install(FILES ${python_files} DESTINATION lib/arv)
