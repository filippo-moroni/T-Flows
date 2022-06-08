!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   Called at the beginning of simulation and it has four distinct parts:      !
!     1. Read physical properties for air                                      !
!     2. Impose initial temperature field                                      !
!     3. Set initial physical properties                                       !
!     4. Compute and store buoyancy forces due to this density variation       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! current time step
  real,    intent(in)         :: time     ! current physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: t
  integer                  :: i, c, fu
  real                     :: wi, wip
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Energy(t)

  !-----------------------------------------!
  !                                         !
  !   1. Read physical properties for air   !
  !                                         !
  !-----------------------------------------!

  !----------------------------------------!
  !   Open file with physical properties   !
  !----------------------------------------!
  call File % Open_For_Reading_Ascii("air_properties_at_1_bar.dat", fu)

  !-----------------------------!
  !   Read all the properties   !
  !-----------------------------!
  do i = 1, N_ITEMS
    call File % Read_Line(fu)

    ! Read desired properties
    read(line % tokens(1), *) air_t(i)
    read(line % tokens(2), *) air_rho(i)
    read(line % tokens(3), *) air_mu(i)
    read(line % tokens(5), *) air_cp(i)
    read(line % tokens(6), *) air_lambda(i)

    ! Fix units where needed (check the values in the table)
    air_cp(i) = air_cp(i) * 1.0e3
    air_mu(i) = air_mu(i) / 1.0e5
  end do

  close(fu)

  if(this_proc < 2) then
    print '(a)',        ' #============================================'
    print '(a)',        ' # Output from user function, read properties!'
    print '(a)',        ' #--------------------------------------------'
  end if

  !-----------------------------------!
  !                                   !
  !   2. Impose initial temperature   !
  !                                   !
  !-----------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    t % n(c)  = -2.0 + 0.008 * Grid % wall_dist(c)
    t % o(c)  = t % n(c)
    t % oo(c) = t % n(c)
  end do

  !----------------------------------------!
  !                                        !
  !   3. Set initial physical properties   !
  !                                        !
  !----------------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells

    ! Browse through all table entries
    do i = 1, N_ITEMS - 1

      ! Did you find the right interval
      if(Flow % t % n(c) >= air_t(i) .and.  &
         Flow % t % n(c) <  air_t(i+1)) then

        ! If so, calculate interpolation factors ...
        wi  = (air_t(i+1) - Flow % t % n(c)) / (air_t(i+1) - air_t(i))
        wip = 1.0 - wi

        ! ... and interpolate physical properties
        Flow % density(c)      = wi * air_rho   (i)  + wip * air_rho   (i+1)
        Flow % viscosity(c)    = wi * air_mu    (i)  + wip * air_mu    (i+1)
        Flow % conductivity(c) = wi * air_lambda(i)  + wip * air_lambda(i+1)
        Flow % capacity(c)     = wi * air_cp    (i)  + wip * air_cp    (i+1)
      end if

    end do
  end do

  !------------------------------------------------------------------------!
  !                                                                        !
  !   4. Compute and store buoyancy forces due to this density variation   !
  !                                                                        !
  !------------------------------------------------------------------------!
  allocate(static_cell_fx(-Grid % n_bnd_cells:Grid % n_cells))
  allocate(static_cell_fy(-Grid % n_bnd_cells:Grid % n_cells))
  allocate(static_cell_fz(-Grid % n_bnd_cells:Grid % n_cells))
  call Flow % Buoyancy_Forces(1);  static_cell_fx(:) = Flow % cell_fx(:)
  call Flow % Buoyancy_Forces(2);  static_cell_fy(:) = Flow % cell_fy(:)
  call Flow % Buoyancy_Forces(3);  static_cell_fz(:) = Flow % cell_fz(:)

  ! Check what you got
  call Grid % Save_Debug_Vtu("static-forces",                    &
                            vector_cell=(/static_cell_fx,        &
                                          static_cell_fy,        &
                                          static_cell_fz/),      &
                            vector_name="Static Force [N/m^3]")

  end subroutine
