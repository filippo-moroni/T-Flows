!==============================================================================!
  subroutine Allocate_Field(Flow, grid)
!------------------------------------------------------------------------------!
!   Allocates memory for the entire field.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type)         :: Flow
  type(Grid_Type),   target :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: sc, nb, nc, nn, nf
  character(VL) :: c_name, q_name
!==============================================================================!

  ! Store the pointer to a grid
  Flow % pnt_grid => grid

  ! Take some aliases
  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes
  nf = grid % n_faces

  !---------------------------------------------!
  !   Allocate memory for physical properties   !
  !---------------------------------------------!
  allocate(Flow % density     (-nb:nc));  Flow % density(:)      = 0.0
  allocate(Flow % viscosity   (-nb:nc));  Flow % viscosity(:)    = 0.0
  allocate(Flow % capacity    (-nb:nc));  Flow % capacity(:)     = 0.0
  allocate(Flow % conductivity(-nb:nc));  Flow % conductivity(:) = 0.0

  !----------------------------------!
  !   Memory for gradient matrices   !  (are the latter two used at all?)
  !----------------------------------!
  allocate(Flow % grad_c2c(6, nc));  Flow % grad_c2c(:,:) = 0.0
  allocate(Flow % grad_f2c(6, nc));  Flow % grad_f2c(:,:) = 0.0
  allocate(Flow % grad_n2c(6, nc));  Flow % grad_n2c(:,:) = 0.0
  allocate(Flow % grad_c2n(6, nn));  Flow % grad_c2n(:,:) = 0.0

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Allocate_Solution(Flow % u,   grid, 'U', '')
  call Var_Mod_Allocate_Solution(Flow % v,   grid, 'V', '')
  call Var_Mod_Allocate_Solution(Flow % w,   grid, 'W', '')

  ! Potential for initial velocity computation
  call Var_Mod_Allocate_Solution(Flow % pot, grid, 'POT', '')

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only(Flow % pp, grid, 'PP')
  call Var_Mod_Allocate_New_Only(Flow % p,  grid, 'P')

  ! Allocate memory for volumetric fluxes
  call Face_Mod_Allocate(Flow % v_flux, grid, 'V_FL')

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(Flow % heat_transfer) then
    call Var_Mod_Allocate_Solution(Flow % t, grid, 'T', 'Q')
  end if ! heat_transfer

  allocate(Flow % vort (-nb:nc)); Flow % vort  = 0.
  allocate(Flow % shear(-nb:nc)); Flow % shear = 0.

  !--------------------------------------------------------------!
  !   Nine variables which follow are needed for Rhie and Chow   !
  !--------------------------------------------------------------!
  allocate(Flow % fx(nc));           Flow % fx = 0.0
  allocate(Flow % fy(nc));           Flow % fy = 0.0
  allocate(Flow % fz(nc));           Flow % fz = 0.0

  allocate(Flow % cell_fx(-nb:nc));  Flow % cell_fx = 0.0
  allocate(Flow % cell_fy(-nb:nc));  Flow % cell_fy = 0.0
  allocate(Flow % cell_fz(-nb:nc));  Flow % cell_fz = 0.0

  allocate(Flow % face_fx(nf));      Flow % face_fx = 0.0
  allocate(Flow % face_fy(nf));      Flow % face_fy = 0.0
  allocate(Flow % face_fz(nf));      Flow % face_fz = 0.0

  !--------------------------------------!
  !   Allocate memory for user scalars   !
  !--------------------------------------!
  allocate(Flow % scalar(Flow % n_scalars))

  !-------------------------------------!
  !   Browse through all user scalars   !
  !-------------------------------------!
  do sc = 1, Flow % n_scalars

    ! Set variable name
    c_name = 'C_00'
    q_name = 'Q_00'
    write(c_name(3:4),'(i2.2)') sc
    write(q_name(3:4),'(i2.2)') sc

    ! Allocate memory for passive scalar
    call Var_Mod_Allocate_Solution(Flow % scalar(sc), grid, c_name, q_name)

  end do

  end subroutine