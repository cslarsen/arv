cmake_minimum_required(VERSION 3.5)
project(arv-tests)

enable_testing()

add_test(arv-python
  ${PYTHON_EXECUTABLE} -m unittest discover ${CMAKE_SOURCE_DIR}/tests -v)

add_test(arv-benchmark
  ${PYTHON_EXECUTABLE}
  ${CMAKE_SOURCE_DIR}/tests/bench.py
  ${CMAKE_SOURCE_DIR}/tests/fake_genome.txt -c10)

set(env
  "PYTHONPATH=${CMAKE_BINARY_DIR}/cython:${CMAKE_BINARY_DIR}/python"
  "GENOME=${CMAKE_SOURCE_DIR}/tests/fake_genome.txt")

set_tests_properties(arv-python PROPERTIES ENVIRONMENT "${env}")
set_tests_properties(arv-benchmark PROPERTIES ENVIRONMENT "${env}")
