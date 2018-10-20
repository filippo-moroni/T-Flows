!==============================================================================!
  real function Y_Plus_Low_Re(u_tau, wall_dist, kin_vis)
!------------------------------------------------------------------------------!
!   Calculates y+ for low Reynolds approach.                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: TINY
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: u_tau, wall_dist, kin_vis
!==============================================================================!

  Y_Plus_Low_Re = max(u_tau * wall_dist / kin_vis, TINY)

  end function
