include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  real,            pointer :: dt
  integer                  :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  dt   => Flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  fun % n(:) = 0.0

  ! Ellipsoid:
  call Vof_Initialization_Ellipsoid(Vof)

  ! Naive way to update bounary values
  do s = 1, Grid % n_bnd_cells
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      fun % n(c2) = fun % n(c1)
    end if
  end do

  ! Update buffer values
  call Grid % Exchange_Cells_Real(fun % n)

  ! Set old values to be the same as new ones
  fun % o(:) = fun % n(:)

  !--------------------!
  !   Create a front   !
  !--------------------!
  if(Vof % track_front) then
    call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                            .false.)  ! don't print messages
    call Vof % Front % Print_Front_Statistics()
  end if

  end subroutine