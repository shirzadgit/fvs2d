subroutine initialize(t)

  use data_grid,    only  : ncells
  use data_sol
  use mms_euler2d,  only  : mms_sol

  implicit none
  real,intent(in)   :: t


  !----------------------------------------------------------------------------!
  ! initialize
  !----------------------------------------------------------------------------!
  pvar(1:nvar,1:ncells) = 0.d0


  !----------------------------------------------------------------------------!
  ! initialize with manufactured solution
  !----------------------------------------------------------------------------!
  pvar(1:nvar,1:ncells) = mms_sol(1:nvar,1:ncells)


  return
end subroutine initialize
