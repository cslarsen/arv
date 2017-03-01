cmake_minimum_required(VERSION 3.5)
project(arv-tests)

enable_testing()

add_test(testarvcy.py COMMAND
  ${PYTHON_EXECUTABLE} -m unittest discover ${CMAKE_CURRENT_BINARY_DIR} -v)

configure_file(fake_genome.txt fake_genome.txt COPYONLY)
configure_file(testarvcy.py testarvcy.py COPYONLY)

set_tests_properties(testarvcy.py PROPERTIES
  ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../cython")
