!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                             curr_dt, time)  ! name shouldn't have a minus sign
!------------------------------------------------------------------------------!
!   Calculates the non-linear, tensorial "viscosity", and the associated SGS   !
!   stress tensor.                                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c
  real                      :: ixp, iyp, izp, ixyp, ixzp, iyzp
  ! Variable names should be all in small leters,
  ! unless they are classes like Flow, Grid, Turb ...
  real                      :: du_dx, du_dy, du_dz,  &   ! line should be max 80
                               dv_dx, dv_dy, dv_dz,  &
                               dw_dx, dw_dy, dw_dz
  real                      :: sgv_x, sgv_y, sgv_z,     &   ! line should be max 80
                               sgv_xy, sgv_yx, sgv_xz,  &
                               sgv_zx, sgv_yz, sgv_zy
!------------------------------------------------------------------------------!
!                                                                              !
!   nii_ki = (1/2*Vol) * I'_kh * dUi/dxh                                       !
!                                                                              !
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Initialize the tensors
  nlin_sgv % x (:) = 0.0
  nlin_sgv % y (:) = 0.0
  nlin_sgv % z (:) = 0.0
  nlin_sgv % xy(:) = 0.0
  nlin_sgv % yx(:) = 0.0
  nlin_sgv % xz(:) = 0.0
  nlin_sgv % zx(:) = 0.0
  nlin_sgv % yz(:) = 0.0
  nlin_sgv % zy(:) = 0.0

  nlin_stress % x (:) = 0.0
  nlin_stress % y (:) = 0.0
  nlin_stress % z (:) = 0.0
  nlin_stress % xy(:) = 0.0
  nlin_stress % yx(:) = 0.0
  nlin_stress % xz(:) = 0.0
  nlin_stress % zx(:) = 0.0
  nlin_stress % yz(:) = 0.0
  nlin_stress % zy(:) = 0.0
  
  do c = 1, Grid % n_cells
    ! Take aliases
    ixp = ( 0.5 * (  Grid % ixx (c)  &
                   - Grid % iyy (c)  &
                   - Grid % izz (c) ) )     ! I'xx         (I'11)
    iyp = ( 0.5 * (- Grid % ixx (c)  &
                   + Grid % iyy (c)  &
                   - Grid % izz (c) ) )     ! I'yy         (I'22)
    izp = ( 0.5 * (- Grid % ixx (c)  &
                   - Grid % iyy (c)  &
                   + Grid % izz (c) ) )     ! I'zz         (I'33)
    ixyp = Grid % ixy (c)              ! I'xy = I'yx  (I'12 = I'21)
    ixzp = Grid % ixz (c)              ! I'xz = I'zx  (I'13 = I'31)
    iyzp = Grid % iyz (c)              ! I'yz = I'zy  (I'23 = I'32)

    du_dx = u % x(c)  ! dU/dx    ¿Podría estar aquí el error? Por estar cogiendo el gradiente de un campo de velocidad que no sea el correcto (presumiblemente)
    du_dy = u % y(c)  ! dU/dy
    du_dz = u % z(c)  ! dU/dz
    dv_dx = v % x(c)  ! dV/dx
    dv_dy = v % y(c)  ! dV/dy
    dv_dz = v % z(c)  ! dV/dz
    dw_dx = w % x(c)  ! dW/dx
    dw_dy = w % y(c)  ! dW/dy
    dw_dz = w % z(c)  ! dW/dz

    ! Construct the non-linear "viscosity" tensor
    nlin_sgv % x  (c) = (0.5 / Grid % vol(c)) * ( ixp  * du_dx  &
                                                + ixyp * du_dy  &
                                                + ixzp * du_dz )
    nlin_sgv % y  (c) = (0.5 / Grid % vol(c)) * ( ixyp * dv_dx  &
                                                + iyp  * dv_dy  &
                                                + iyzp * dv_dz )
    nlin_sgv % z  (c) = (0.5 / Grid % vol(c)) * ( ixzp * dw_dx  &
                                                + iyzp * dw_dy  &
                                                + izp  * dw_dz )
    nlin_sgv % xy (c) = (0.5 / Grid % vol(c)) * ( ixp  * dv_dx  &
                                                + ixyp * dv_dy  &
                                                + ixzp * dv_dz )
    nlin_sgv % yx (c) = (0.5 / Grid % vol(c)) * ( ixyp * du_dx  &
                                                + iyp  * du_dy  &
                                                + iyzp * du_dz )
    nlin_sgv % xz (c) = (0.5 / Grid % vol(c)) * ( ixp  * dw_dx  &
                                                + ixyp * dw_dy  &
                                                + ixzp * dw_dz )
    nlin_sgv % zx (c) = (0.5 / Grid % vol(c)) * ( ixzp * du_dx  &
                                                + iyzp * du_dy  &
                                                + izp  * du_dz )
    nlin_sgv % yz (c) = (0.5 / Grid % vol(c)) * ( ixyp * dw_dx  &
                                                + iyp  * dw_dy  &
                                                + iyzp * dw_dz )
    nlin_sgv % zy (c) = (0.5 / Grid % vol(c)) * ( ixzp * dv_dx  &
                                                + iyzp * dv_dy  &
                                                + izp  * dv_dz )

    ! Since we're at it, we may as well construct the subgrid STRESS tensor too
      ! Take aliases
      sgv_x  = nlin_sgv % x  (c)
      sgv_y  = nlin_sgv % y  (c)
      sgv_z  = nlin_sgv % z  (c)
      sgv_xy = nlin_sgv % xy (c)
      sgv_yx = nlin_sgv % yx (c)
      sgv_xz = nlin_sgv % xz (c)
      sgv_zx = nlin_sgv % zx (c)
      sgv_yz = nlin_sgv % yz (c)
      sgv_zy = nlin_sgv % zy (c)

      ! Now actually build the thing
      nlin_stress % x  (c) = (- sgv_x  * du_dx - sgv_x  * du_dx)  &
                           + (- sgv_yx * du_dy - sgv_yx * du_dy)  &
                           + (- sgv_zx * du_dz - sgv_zx * du_dz)

      nlin_stress % y  (c) = (- sgv_xy * dv_dx - sgv_xy * dv_dx)  &
                           + (- sgv_y  * dv_dy - sgv_y  * dv_dy)  &
                           + (- sgv_zy * dv_dz - sgv_zy * dv_dz)

      nlin_stress % z  (c) = (- sgv_xz * dw_dx - sgv_xz * dw_dx)  &
                           + (- sgv_yz * dw_dy - sgv_yz * dw_dy)  &
                           + (- sgv_z  * dw_dz - sgv_z  * dw_dz)

      nlin_stress % xy (c) = (- sgv_xy * du_dx - sgv_x  * dv_dx)  &
                           + (- sgv_y  * du_dy - sgv_yx * dv_dy)  &
                           + (- sgv_zy * du_dz - sgv_zx * dv_dz)

      nlin_stress % yx (c) = (- sgv_x  * dv_dx - sgv_xy * du_dx)  &
                           + (- sgv_yx * dv_dy - sgv_y  * du_dy)  &
                           + (- sgv_zx * dv_dz - sgv_zy * du_dz)

      nlin_stress % xz (c) = (- sgv_xz * du_dx - sgv_x  * dw_dx)  &
                           + (- sgv_yz * du_dy - sgv_yx * dw_dy)  &
                           + (- sgv_z  * du_dz - sgv_zx * dw_dz)

      nlin_stress % zx (c) = (- sgv_x  * dw_dx - sgv_xz * du_dx)  &
                           + (- sgv_yx * dw_dy - sgv_yz * du_dy)  &
                           + (- sgv_zx * dw_dz - sgv_z  * du_dz)

      nlin_stress % yz (c) = (- sgv_xz * dv_dx - sgv_xy * dw_dx)  &
                           + (- sgv_yz * dv_dy - sgv_y  * dw_dy)  &
                           + (- sgv_z  * dv_dz - sgv_zy * dw_dz)

      nlin_stress % zy (c) = (- sgv_xy * dw_dx - sgv_xz * dv_dx)  &
                           + (- sgv_y  * dw_dy - sgv_yz * dv_dy)  &
                           + (- sgv_zy * dw_dz - sgv_z  * dv_dz)

  end do
  
  
  print *, 'SGS non-linear viscosity and stress tensors calculated'
                             
  call Grid % Save_Debug_Vtu('niiSGS_ij',                      &
                             tensor_cell = (/nlin_sgv % x,     &
                                             nlin_sgv % y,     &
                                             nlin_sgv % z,     &
                                             nlin_sgv % xy,    &
                                             nlin_sgv % yx,    &
                                             nlin_sgv % xz,    &
                                             nlin_sgv % zx,    &
                                             nlin_sgv % yz,    &
                                             nlin_sgv % zy/),  &
                             tensor_comp = 9,                  &
                             tensor_name = 'Nonlinear Turbulent Viscosity Tensor')
                             
  call Grid % Save_Debug_Vtu('tauSGS_ij',                         &
                             tensor_cell = (/nlin_stress % x,     &
                                             nlin_stress % y,     &
                                             nlin_stress % z,     &
                                             nlin_stress % xy,    &
                                             nlin_stress % yx,    &
                                             nlin_stress % xz,    &
                                             nlin_stress % zx,    &
                                             nlin_stress % yz,    &
                                             nlin_stress % zy/),  &
                             tensor_comp = 9,                     &
                             tensor_name = 'Nonlinear Turbulent Viscosity Tensor')

  end subroutine
