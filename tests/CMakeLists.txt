cmake_minimum_required(VERSION 3.5)
project(arv-tests)

enable_testing()

file(GLOB python_files
  "${CMAKE_SOURCE_DIR}/tests/*.py")

foreach(filename ${python_files})
  get_filename_component(basename ${filename} NAME)
  configure_file(${filename} ${CMAKE_BINARY_DIR}/tests/${basename} COPYONLY)
endforeach()

add_test(arv-python
  ${PYTHON_EXECUTABLE} -m unittest discover ${CMAKE_BINARY_DIR}/tests -v)

add_test(arv-benchmark
  ${PYTHON_EXECUTABLE}
  ${CMAKE_BINARY_DIR}/tests/bench.py
  ${CMAKE_SOURCE_DIR}/tests/fake_genome.txt -c10)

set(env
  "PYTHONPATH=${CMAKE_BINARY_DIR}/cython:${CMAKE_BINARY_DIR}/python"
  "GENOME=${CMAKE_SOURCE_DIR}/tests/fake_genome.txt")

set_tests_properties(arv-python PROPERTIES ENVIRONMENT "${env}")
set_tests_properties(arv-benchmark PROPERTIES ENVIRONMENT "${env}")
