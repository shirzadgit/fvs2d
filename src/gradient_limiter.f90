module gradient_limiter

  use mainparam,      only  : nvar
  use input,          only  : lgrad_limiter, lface_reconst_upwind1st
  use data_solution,  only  : grad, phi_lim

  implicit none


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_gradient_limiter

    implicit none


    !--------------------------------------------------------------------------!
    ! if 1storder upwind --> phi_limiter = 0
    !--------------------------------------------------------------------------!
    if (lface_reconst_upwind1st) then
      phi_lim = 0.d0
      return
    endif


    !-- If no limiter
    if (.not.lgrad_limiter) then
      phi_lim = 1.d0
      return
    endif

    !--------------------------------------------------------------------------!
    ! setup
    !--------------------------------------------------------------------------!

    return
  end subroutine compute_gradient_limiter

end module gradient_limiter
