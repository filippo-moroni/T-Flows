!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
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
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Make sure you compiled it
  if(this_proc < 2) print *, '# In User_Mod_Beginning_Of_Simulation'

  ! Take alias
  Grid => Flow % pnt_grid
  
  ! Allocate memory for user arrays
  allocate(nlin_sgv % x (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % y (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % z (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % xy(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % yx(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % xz(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % zx(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % yz(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_sgv % zy(-Grid % n_bnd_cells : Grid % n_cells))

  allocate(nlin_stress % x (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % y (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % z (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % xy(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % yx(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % xz(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % zx(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % yz(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(nlin_stress % zy(-Grid % n_bnd_cells : Grid % n_cells))

  !--------------------------------------------------!
  !   Check the computed cells' moments of inertia   !
  !--------------------------------------------------!
  call Grid % Save_Debug_Vtu('cell-inertia',                &
                             tensor_cell = (/Grid % ixx,    &
                                             Grid % iyy,    &
                                             Grid % izz,    &
                                             Grid % ixy,    &
                                             Grid % ixz,    &
                                             Grid % iyz/),  &
                             tensor_comp = 6,               &
                             tensor_name = 'Cell Inertia')

  end subroutine
