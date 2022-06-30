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
  ! pvar(:)%r=0.d0
  ! pvar(:)%u=0.d0
  ! pvar(:)%v=0.d0
  ! pvar(:)%p=0.d0


  !----------------------------------------------------------------------------!
  ! initialize with manufactured solution
  !----------------------------------------------------------------------------!
  pvar(1:nvar,1:ncells) = mms_sol(1:nvar,1:ncells)
  ! pvar(1:ncells)%r = mms_sol(1:ncells)%r
  ! pvar(1:ncells)%u = mms_sol(1:ncells)%u
  ! pvar(1:ncells)%v = mms_sol(1:ncells)%v
  ! pvar(1:ncells)%p = mms_sol(1:ncells)%p


  return
end subroutine initialize
