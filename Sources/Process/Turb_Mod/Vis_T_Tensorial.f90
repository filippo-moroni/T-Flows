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
!------------------------------------------------------------------------------!
!                                                                              !
! tau_ij = - I'_kh/omega*[(du_j/dx_h)(du_i/dx_k) + (du_i/dx_h)(du_j/dx_k)      !
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

    du_dx = u % x(c)  ! dU/dx
    du_dy = u % y(c)  ! dU/dy
    du_dz = u % z(c)  ! dU/dz
    dv_dx = v % x(c)  ! dV/dx
    dv_dy = v % y(c)  ! dV/dy
    dv_dz = v % z(c)  ! dV/dz
    dw_dx = w % x(c)  ! dW/dx
    dw_dy = w % y(c)  ! dW/dy
    dw_dz = w % z(c)  ! dW/dz

    ! Now actually build the thing
    Turb % tau_11 (c)    = (- sgv_x  * du_dx - sgv_x  * du_dx)  &
                         + (- sgv_yx * du_dy - sgv_yx * du_dy)  &
                         + (- sgv_zx * du_dz - sgv_zx * du_dz)

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

    Turb % tau_13 (c)    = (- sgv_xz * du_dx - sgv_x  * dw_dx)  &
                         + (- sgv_yz * du_dy - sgv_yx * dw_dy)  &
                         + (- sgv_z  * du_dz - sgv_zx * dw_dz)

    Turb % tau_31 (c)    = Turb % tau_13 (c)

    Turb % tau_23 (c)    = (- sgv_xz * dv_dx - sgv_xy * dw_dx)  &
                         + (- sgv_yz * dv_dy - sgv_y  * dw_dy)  &
                         + (- sgv_z  * dv_dz - sgv_zy * dw_dz)

    Turb % tau_32 (c)    = Turb % tau_23 (c)

  end do

  end subroutine
