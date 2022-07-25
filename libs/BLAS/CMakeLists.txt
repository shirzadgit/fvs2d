file(GLOB_RECURSE blas_source_list "*.f")
add_library(blas ${blas_source_list})
target_include_directories (blas PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
