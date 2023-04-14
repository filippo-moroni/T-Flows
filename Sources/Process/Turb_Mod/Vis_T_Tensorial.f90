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
  real                      :: t11, t22, t33
  real                      :: t12, t13, t23
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

  ! Initialize subgrid stress tensor
  Turb % tau_11 = 0.0
  Turb % tau_22 = 0.0
  Turb % tau_33 = 0.0
  Turb % tau_12 = 0.0
  Turb % tau_13 = 0.0
  Turb % tau_23 = 0.0

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
    t11 =     Grid % ixxp(c)*du_dx*du_dx &
        +     Grid % iyyp(c)*du_dy*du_dy &
        +     Grid % izzp(c)*du_dz*du_dz &
        + 2.0*Grid % ixyp(c)*du_dx*du_dy &
        + 2.0*Grid % ixzp(c)*du_dx*du_dz &
        + 2.0*Grid % iyzp(c)*du_dy*du_dz
                           
    t22 =     Grid % ixxp(c)*dv_dx*dv_dx &
        +     Grid % iyyp(c)*dv_dy*dv_dy &
        +     Grid % izzp(c)*dv_dz*dv_dz &
        + 2.0*Grid % ixyp(c)*dv_dx*dv_dy & 
        + 2.0*Grid % ixzp(c)*dv_dx*dv_dz &
        + 2.0*Grid % iyzp(c)*dv_dy*dv_dz
                        
    t33 =     Grid % ixxp(c)*dw_dx*dw_dx &
        +     Grid % iyyp(c)*dw_dy*dw_dy &
        +     Grid % izzp(c)*dw_dz*dw_dz &
        + 2.0*Grid % ixyp(c)*dw_dx*dw_dy &
        + 2.0*Grid % ixzp(c)*dw_dx*dw_dz &
        + 2.0*Grid % iyzp(c)*dw_dy*dw_dz
        
    Turb % tau_11(c) = -1.0/omega*t11
    Turb % tau_22(c) = -1.0/omega*t22
    Turb % tau_33(c) = -1.0/omega*t33
                      
    ! Deviatoric components
    Turb % tau_12(c)  =     Grid % ixxp(c)*(2.0*dv_dx*du_dx)           &
                      +     Grid % iyyp(c)*(2.0*dv_dy*du_dy)           &
                      +     Grid % izzp(c)*(2.0*dv_dz*du_dz)           &
                      + 2.0*Grid % ixyp(c)*(dv_dx*du_dy + dv_dy*du_dx) &
                      + 2.0*Grid % ixzp(c)*(dv_dx*du_dz + dv_dz*du_dx) &
                      + 2.0*Grid % iyzp(c)*(dv_dy*du_dz + dv_dz*du_dy)
                        
    Turb % tau_13(c)  =     Grid % ixxp(c)*(2.0*dw_dx*du_dx)           &
                      +     Grid % iyyp(c)*(2.0*dw_dy*du_dy)           &
                      +     Grid % izzp(c)*(2.0*dw_dz*du_dz)           &
                      + 2.0*Grid % ixyp(c)*(dw_dx*du_dy + dw_dy*du_dx) &
                      + 2.0*Grid % ixzp(c)*(dw_dx*du_dz + dw_dz*du_dx) &
                      + 2.0*Grid % iyzp(c)*(dw_dy*du_dz + dw_dz*du_dy)
                        
    Turb % tau_23(c)  =     Grid % ixxp(c)*(2.0*dw_dx*dv_dx)           &
                      +     Grid % iyyp(c)*(2.0*dw_dy*dv_dy)           &
                      +     Grid % izzp(c)*(2.0*dw_dz*dv_dz)           &
                      + 2.0*Grid % ixyp(c)*(dw_dx*dv_dy + dw_dy*dv_dx) &
                      + 2.0*Grid % ixzp(c)*(dw_dx*dv_dz + dw_dz*dv_dx) &
                      + 2.0*Grid % iyzp(c)*(dw_dy*dv_dz + dw_dz*dv_dy)
                        
                        
 
                      
                        
                        
                        
    
                        
                        
    

    Turb % tau_22 (c)    = (- sgv_xy * dv_dx - sgv_xy * dv_dx)  &
                         + (- sgv_y  * dv_dy - sgv_y  * dv_dy)  &
                         + (- sgv_zy * dv_dz - sgv_zy * dv_dz)

    Turb % tau_33 (c)    = (- sgv_xz * dw_dx - sgv_xz * dw_dx)  &
                         + (- sgv_yz * dw_dy - sgv_yz * dw_dy)  &
                         + (- sgv_z  * dw_dz - sgv_z  * dw_dz)

    Turb % tau_12 (c)    = (- sgv_xy * du_dx - sgv_x  * dv_dx)  &
                         + (- sgv_y  * du_dy - sgv_yx * dv_dy)  &
                         + (- sgv_zy * du_dz - sgv_zx * dv_dz)

    Turb % tau_21 (c)    = Turb % tau_12 (c)

    Turb % tau_31 (c)    = Turb % tau_13 (c)

    Turb % tau_32 (c)    = Turb % tau_23 (c)

  end do

  end subroutine
