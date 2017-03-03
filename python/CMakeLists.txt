cmake_minimum_required(VERSION 3.5)
project(arv-python CXX)

configure_file(setup.py.in setup.py)

file(GLOB python_files
  ${CMAKE_CURRENT_SOURCE_DIR}/arv/*.py
  ${CMAKE_CURRENT_BUILD_DIR}/setup.py
)

# Copy to build folder so we can run tests
foreach(filename ${python_files})
  get_filename_component(basename ${filename} NAME)
  configure_file(${filename} arv/${basename} COPYONLY)
endforeach()

install(FILES ${python_files} DESTINATION lib/arv)
