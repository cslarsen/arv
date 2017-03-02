cmake_minimum_required(VERSION 3.5)
project(sparsehash-config)

# Run configure to produce config.h, which will be copied as
# cpp/sparsehash/internal/sparseconfig.h
add_custom_command(
  OUTPUT "${PROJECT_SOURCE_DIR}/sparsehash/src/config.h"
  COMMAND "${PROJECT_SOURCE_DIR}/sparsehash/configure"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/sparsehash"
  USES_TERMINAL
  COMMENT "Running sparsehash/configure"
)
