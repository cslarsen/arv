include(FindPackageHandleStandardArgs)

find_path(google_densehash_INCLUDE_DIR
          google/dense_hash_map)

find_package_handle_standard_args(
  google_densehash DEFAULT_MSG
  google_densehash_INCLUDE_DIR)
