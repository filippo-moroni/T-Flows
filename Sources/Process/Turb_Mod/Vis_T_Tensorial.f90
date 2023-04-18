!==============================================================================!
  subroutine Vis_T_Tensorial(Turb)
!------------------------------------------------------------------------------!
! Calculates the subgrid stress tensor of TVM model.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c
  real                      :: du_dx, du_dy, du_dz
  real                      :: dv_dx, dv_dy, dv_dz
  real                      :: dw_dx, dw_dy, dw_dz
  real                      :: omega
!------------------------------------------------------------------------------!
!                                                                              !
!   tau_ij = - I'_kh/omega*[(du_j/dx_h)(du_i/dx_k) + (du_i/dx_h)(du_j/dx_k)    !
!                                                                              !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  
  ! Initialize
  Turb % tau_11 = 0.0
  Turb % tau_22 = 0.0
  Turb % tau_33 = 0.0
  Turb % tau_12 = 0.0
  Turb % tau_13 = 0.0
  Turb % tau_23 = 0.0
  
  ! Cycle through all cells
  do c = 1, Grid % n_cells
  
    ! Take aliases for velocity components gradients 
    du_dx = u % x(c)  ! dU/dx
    du_dy = u % y(c)  ! dU/dy
    du_dz = u % z(c)  ! dU/dz
    dv_dx = v % x(c)  ! dV/dx
    dv_dy = v % y(c)  ! dV/dy
    dv_dz = v % z(c)  ! dV/dz
    dw_dx = w % x(c)  ! dW/dx
    dw_dy = w % y(c)  ! dW/dy
    dw_dz = w % z(c)  ! dW/dz
    
    ! Take alias for grid element volume
    omega = Grid % vol(c)
    
    ! Diagonal components
    Turb % tau_11(c) = - ( Grid % ixxp(c)*du_dx*du_dx   &
                     +     Grid % iyyp(c)*du_dy*du_dy   &
                     +     Grid % izzp(c)*du_dz*du_dz   &
                     + 2.0*Grid % ixyp(c)*du_dx*du_dy   &
                     + 2.0*Grid % ixzp(c)*du_dx*du_dz   &
                     + 2.0*Grid % iyzp(c)*du_dy*du_dz)) &
                     /omega
                           
    Turb % tau_22(c) =   -(Grid % ixxp(c)*dv_dx*dv_dx   &
                     +     Grid % iyyp(c)*dv_dy*dv_dy   &
                     +     Grid % izzp(c)*dv_dz*dv_dz   &
                     + 2.0*Grid % ixyp(c)*dv_dx*dv_dy   & 
                     + 2.0*Grid % ixzp(c)*dv_dx*dv_dz   &
                     + 2.0*Grid % iyzp(c)*dv_dy*dv_dz)) &
                     /omega
                        
    Turb % tau_33(c) =   -(Grid % ixxp(c)*dw_dx*dw_dx   &
                     +     Grid % iyyp(c)*dw_dy*dw_dy   &
                     +     Grid % izzp(c)*dw_dz*dw_dz   &
                     + 2.0*Grid % ixyp(c)*dw_dx*dw_dy   &
                     + 2.0*Grid % ixzp(c)*dw_dx*dw_dz   &
                     + 2.0*Grid % iyzp(c)*dw_dy*dw_dz)) &
                     /omega
                         
    ! Deviatoric components
    Turb % tau_12(c) = -(Grid % ixxp(c)*(dv_dx*du_dx)                &
                       + Grid % iyyp(c)*(dv_dy*du_dy)                &
                       + Grid % izzp(c)*(dv_dz*du_dz)                &
                       + Grid % ixyp(c)*(dv_dx*du_dy + dv_dy*du_dx)  &
                       + Grid % ixzp(c)*(dv_dx*du_dz + dv_dz*du_dx)  & 
                       + Grid % iyzp(c)*(dv_dy*du_dz + dv_dz*du_dy)) &
                       /omega
                        
    Turb % tau_13(c) = -(Grid % ixxp(c)*(dw_dx*du_dx)                &
                       + Grid % iyyp(c)*(dw_dy*du_dy)                &
                       + Grid % izzp(c)*(dw_dz*du_dz)                &
                       + Grid % ixyp(c)*(dw_dx*du_dy + dw_dy*du_dx)  &
                       + Grid % ixzp(c)*(dw_dx*du_dz + dw_dz*du_dx)  &
                       + Grid % iyzp(c)*(dw_dy*du_dz + dw_dz*du_dy)) &
                       /omega
                        
    Turb % tau_23(c) = -(Grid % ixxp(c)*(dw_dx*dv_dx)                &
                       + Grid % iyyp(c)*(dw_dy*dv_dy)                &
                       + Grid % izzp(c)*(dw_dz*dv_dz)                &
                       + Grid % ixyp(c)*(dw_dx*dv_dy + dw_dy*dv_dx)  &
                       + Grid % ixzp(c)*(dw_dx*dv_dz + dw_dz*dv_dx)  &
                       + Grid % iyzp(c)*(dw_dy*dv_dz + dw_dz*dv_dy)) &
                       /omega
        
  end do

  end subroutine
