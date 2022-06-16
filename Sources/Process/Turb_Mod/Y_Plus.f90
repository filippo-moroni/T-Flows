!==============================================================================!
  real function Y_Plus(turb, u_tau, wall_dist, kin_vis, z_o)
!------------------------------------------------------------------------------!
!   Calculates y+  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turb_Mod, only: Turb_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
  real            :: u_tau, wall_dist, kin_vis, z_o
!==============================================================================!

  Y_Plus = u_tau * (wall_dist + z_o) / kin_vis

  end function
