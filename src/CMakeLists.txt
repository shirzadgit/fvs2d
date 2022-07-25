include_directories(${CMAKE_BINARY_DIR}/src/)

add_executable(fvs2d_nlin.exe
  fvs2d.f90
  mainparam.f90
  input.f90
  kdtree2.f90
  ios_unstrc.f90
  data_grid.f90
  data_solution.f90
  grid_procs.f90
  interpolation.f90
  gradient_ggcb.f90
  gradient_ggnb.f90
  gradient_lsq.f90
  gradient.f90
  gradient_limiter.f90
  mms.f90
  flux_invscid.f90
  residual.f90
  runge_kutta.f90
  io.f90
  initialize.f90
  test.f90)
          
target_link_libraries (fvs2d_nlin.exe LINK_PUBLIC blas)
