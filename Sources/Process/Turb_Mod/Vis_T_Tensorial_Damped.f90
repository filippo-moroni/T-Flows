!==============================================================================!
  subroutine Vis_T_Tensorial_Damped(Turb)
!------------------------------------------------------------------------------!
!   Calculates the non-linear, tensorial "viscosity", and the associated SGS   !
!   stress tensor, with Van Driest damping at the wall.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: A_POW = 8.3
  integer, parameter :: B_POW = 1.0/7.0
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c
  real                      :: ixp, iyp, izp, ixyp, ixzp, iyzp
  real                      :: du_dx, du_dy, du_dz
  real                      :: dv_dx, dv_dy, dv_dz
  real                      :: dw_dx, dw_dy, dw_dz
  real                      :: sgv_x , sgv_y , sgv_z
  real                      :: sgv_xy, sgv_yx, sgv_xz
  real                      :: sgv_zx, sgv_yz, sgv_zy
  real                      :: nu, dely, u_tan, u_tau, dc
!------------------------------------------------------------------------------!
!                                                                              !
!   nii_ki = (1/2*Vol) * I'_kh * dUi/dxh                                       !
!                                                                              !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Initialize the tensors
  Turb % ten_turb_11 = 0.0
  Turb % ten_turb_12 = 0.0
  Turb % ten_turb_13 = 0.0
  Turb % ten_turb_21 = 0.0
  Turb % ten_turb_22 = 0.0
  Turb % ten_turb_23 = 0.0
  Turb % ten_turb_31 = 0.0
  Turb % ten_turb_32 = 0.0
  Turb % ten_turb_33 = 0.0

  Turb % tau_11 = 0.0
  Turb % tau_12 = 0.0
  Turb % tau_13 = 0.0
  Turb % tau_21 = 0.0
  Turb % tau_22 = 0.0
  Turb % tau_23 = 0.0
  Turb % tau_31 = 0.0
  Turb % tau_32 = 0.0
  Turb % tau_33 = 0.0

  do c = 1, Grid % n_cells
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

    du_dx = u % x(c)  ! dU/dx
    du_dy = u % y(c)  ! dU/dy
    du_dz = u % z(c)  ! dU/dz
    dv_dx = v % x(c)  ! dV/dx
    dv_dy = v % y(c)  ! dV/dy
    dv_dz = v % z(c)  ! dV/dz
    dw_dx = w % x(c)  ! dW/dx
    dw_dy = w % y(c)  ! dW/dy
    dw_dz = w % z(c)  ! dW/dz

    ! Construct the non-linear "viscosity" tensor
    Turb % ten_turb_11 (c) = (0.5 / Grid % vol(c)) * ( ixp  * du_dx  &
                                                     + ixyp * du_dy  &
                                                     + ixzp * du_dz )
    Turb % ten_turb_22 (c) = (0.5 / Grid % vol(c)) * ( ixyp * dv_dx  &
                                                     + iyp  * dv_dy  &
                                                     + iyzp * dv_dz )
    Turb % ten_turb_33 (c) = (0.5 / Grid % vol(c)) * ( ixzp * dw_dx  &
                                                     + iyzp * dw_dy  &
                                                     + izp  * dw_dz )
    Turb % ten_turb_12 (c) = (0.5 / Grid % vol(c)) * ( ixp  * dv_dx  &
                                                     + ixyp * dv_dy  &
                                                     + ixzp * dv_dz )
    Turb % ten_turb_21 (c) = (0.5 / Grid % vol(c)) * ( ixyp * du_dx  &
                                                     + iyp  * du_dy  &
                                                     + iyzp * du_dz )
    Turb % ten_turb_13 (c) = (0.5 / Grid % vol(c)) * ( ixp  * dw_dx  &
                                                     + ixyp * dw_dy  &
                                                     + ixzp * dw_dz )
    Turb % ten_turb_31 (c) = (0.5 / Grid % vol(c)) * ( ixzp * du_dx  &
                                                     + iyzp * du_dy  &
                                                     + izp  * du_dz )
    Turb % ten_turb_23 (c) = (0.5 / Grid % vol(c)) * ( ixyp * dw_dx  &
                                                     + iyp  * dw_dy  &
                                                     + iyzp * dw_dz )
    Turb % ten_turb_32 (c) = (0.5 / Grid % vol(c)) * ( ixzp * dv_dx  &
                                                     + iyzp * dv_dy  &
                                                     + izp  * dv_dz )

    ! Aliases for tensorial viscosity
    sgv_x  = Turb % ten_turb_11 (c)
    sgv_y  = Turb % ten_turb_22 (c)
    sgv_z  = Turb % ten_turb_33 (c)
    sgv_xy = Turb % ten_turb_12 (c)
    sgv_yx = Turb % ten_turb_21 (c)
    sgv_xz = Turb % ten_turb_13 (c)
    sgv_zx = Turb % ten_turb_31 (c)
    sgv_yz = Turb % ten_turb_23 (c)
    sgv_zy = Turb % ten_turb_32 (c)
    
    !------------------------!
    !   Van Driest damping   !
    !------------------------!
    
    ! Kinematic viscosity, cell wall distance and tangential velocity
    nu   = Flow % viscosity(c) / Flow % density(c)
    dely = Grid % wall_dist(c)
    u_tan = sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2)
    
    ! Friction velocity and y+ 
    u_tau = (u_tan/A_POW * (nu/dely)**B_POW) ** (1.0/(1.0+B_POW))
    Turb % y_plus(c) = Grid % wall_dist(c) * u_tau / Flow % viscosity(c)
    
    ! Damping coefficient
    dc = (1.0 - exp(-Turb % y_plus(c) / 25.0))

    !----------------------!
    !   Subgrid stresses   !
    !----------------------!
    
    Turb % tau_11 (c)    = ((- sgv_x  * du_dx - sgv_x  * du_dx)     &
                         +  (- sgv_yx * du_dy - sgv_yx * du_dy)     &
                         +  (- sgv_zx * du_dz - sgv_zx * du_dz))*dc

    Turb % tau_22 (c)    = ((- sgv_xy * dv_dx - sgv_xy * dv_dx)     &
                         +  (- sgv_y  * dv_dy - sgv_y  * dv_dy)     &
                         +  (- sgv_zy * dv_dz - sgv_zy * dv_dz))*dc

    Turb % tau_33 (c)    = ((- sgv_xz * dw_dx - sgv_xz * dw_dx)     &
                         +  (- sgv_yz * dw_dy - sgv_yz * dw_dy)     &
                         +  (- sgv_z  * dw_dz - sgv_z  * dw_dz))*dc

    Turb % tau_12 (c)    = ((- sgv_xy * du_dx - sgv_x  * dv_dx)     &
                         +  (- sgv_y  * du_dy - sgv_yx * dv_dy)     &
                         +  (- sgv_zy * du_dz - sgv_zx * dv_dz))*dc

    Turb % tau_21 (c)    = ((- sgv_x  * dv_dx - sgv_xy * du_dx)     &
                         +  (- sgv_yx * dv_dy - sgv_y  * du_dy)     &
                         +  (- sgv_zx * dv_dz - sgv_zy * du_dz)*dc

    Turb % tau_13 (c)    = ((- sgv_xz * du_dx - sgv_x  * dw_dx)     &
                         +  (- sgv_yz * du_dy - sgv_yx * dw_dy)     &
                         +  (- sgv_z  * du_dz - sgv_zx * dw_dz))*dc

    Turb % tau_31 (c)    = ((- sgv_x  * dw_dx - sgv_xz * du_dx)     &
                         +  (- sgv_yx * dw_dy - sgv_yz * du_dy)     &
                         +  (- sgv_zx * dw_dz - sgv_z  * du_dz))*dc

    Turb % tau_23 (c)    = ((- sgv_xz * dv_dx - sgv_xy * dw_dx)     &
                         +  (- sgv_yz * dv_dy - sgv_y  * dw_dy)     &
                         +  (- sgv_z  * dv_dz - sgv_zy * dw_dz))*dc

    Turb % tau_32 (c)    = ((- sgv_xy * dw_dx - sgv_xz * dv_dx)     &
                         +  (- sgv_y  * dw_dy - sgv_yz * dv_dy)     &
                         +  (- sgv_zy * dw_dz - sgv_z  * dv_dz))*dc

  end do

  end subroutine
