!==============================================================================!
!   Introduce new types to be used with User_Mod                               !
!------------------------------------------------------------------------------!

  ! For physical properties
  integer, parameter :: N_ITEMS = 30

  real :: air_t     (N_ITEMS)
  real :: air_rho   (N_ITEMS)
  real :: air_mu    (N_ITEMS)
  real :: air_cp    (N_ITEMS)
  real :: air_lambda(N_ITEMS)

  ! For hydrostatic force
  real, allocatable :: static_cell_fx(:)
  real, allocatable :: static_cell_fy(:)
  real, allocatable :: static_cell_fz(:)
