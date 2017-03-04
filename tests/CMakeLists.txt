cmake_minimum_required(VERSION 3.5)
project(arv-tests)

enable_testing()

add_test(arv-python
  ${PYTHON_EXECUTABLE} -m unittest discover ${CMAKE_CURRENT_BINARY_DIR} -v)

add_test(arv-benchmark
  ${PYTHON_EXECUTABLE} bench.py fake_genome.txt -c10)

configure_file(bench.py bench.py COPYONLY)
configure_file(fake_genome.txt fake_genome.txt COPYONLY)
configure_file(testarv.py testarv.py COPYONLY)
configure_file(testarvcy.py testarvcy.py COPYONLY)
configure_file(testinfer.py testinfer.py COPYONLY)

set_tests_properties(arv-python PROPERTIES
  ENVIRONMENT
    "PYTHONPATH=${CMAKE_BINARY_DIR}/cython:${CMAKE_BINARY_DIR}/python")

set_tests_properties(arv-benchmark PROPERTIES
  ENVIRONMENT
    "PYTHONPATH=${CMAKE_BINARY_DIR}/cython:${CMAKE_BINARY_DIR}/python")
