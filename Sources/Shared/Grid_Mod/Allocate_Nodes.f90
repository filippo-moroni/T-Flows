!==============================================================================!
  subroutine Grid_Mod_Allocate_Nodes(grid, nn)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nn
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
!==============================================================================!

  ! Store the number of nodes for the grid
  grid % n_nodes = nn

  ! Allocate memory for node coordinates
  allocate(grid % xn(1:nn));  grid % xn(:) = 0.0
  allocate(grid % yn(1:nn));  grid % yn(:) = 0.0
  allocate(grid % zn(1:nn));  grid % zn(:) = 0.0

  allocate(grid % comm % node_glo(1:nn));  grid % comm % node_glo(:) = 0
  do n = 1, nn
    grid % comm % node_glo(n) = n
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  allocate(grid % new_n(n));  grid % new_n(:) = 0
  allocate(grid % old_n(n));  grid % old_n(:) = 0

  end subroutine
