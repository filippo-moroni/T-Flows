!==============================================================================!
  subroutine Vof_Initialization_Plane(Vof)
!------------------------------------------------------------------------------!
!   Initialize as vof = 1 all cells beneath the plane given by 3 points        !
!   sorted anticlockwise to the direction of the plane normal                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, allocatable         :: p_xyz(:,:)
  integer                   :: c, n, fu
  integer                   :: ee, n_planes
  real                      :: n_xyz(3), v1aux(3), v2aux(3)
  real                      :: res_dummy
  real                      :: dd
  real, contiguous, pointer :: prelim_vof(:),min_max_crit_1(:),min_max_crit_2(:)
!==============================================================================!

  call Work % Connect_Real_Cell(prelim_vof, min_max_crit_1, min_max_crit_2)

  ! First take aliases
  Grid => Vof % pnt_grid

  prelim_vof = 0.0

  ! Open file to read Plane parameters:
  call File % Open_For_Reading_Ascii('plane_parameters.ini', fu)

  call File % Read_Line(fu)
  read(line % tokens(1), *) n_planes

  if (allocated(p_xyz)) deallocate(p_xyz)
  allocate(p_xyz(n_planes*3,3))

  do ee = 1, n_planes

    ! Taking extrema
    min_max_crit_1(:) =  HUGE
    min_max_crit_2(:) = -HUGE

    ! Read point 1
    call File % Read_Line(fu)
    read(line % tokens(1), *) p_xyz(1, 1)
    read(line % tokens(2), *) p_xyz(1, 2)
    read(line % tokens(3), *) p_xyz(1, 3)

    ! Read point 2
    call File % Read_Line(fu)
    read(line % tokens(1), *) p_xyz(2, 1)
    read(line % tokens(2), *) p_xyz(2, 2)
    read(line % tokens(3), *) p_xyz(2, 3)

    ! Read point 3
    call File % Read_Line(fu)
    read(line % tokens(1), *) p_xyz(3, 1)
    read(line % tokens(2), *) p_xyz(3, 2)
    read(line % tokens(3), *) p_xyz(3, 3)

    v1aux(:) = p_xyz(2,:) - p_xyz(1,:)
    v2aux(:) = p_xyz(3,:) - p_xyz(1,:)

    n_xyz(1) = v1aux(2) * v2aux(3) - v1aux(3) * v2aux(2)
    n_xyz(2) = - (v1aux(1) * v2aux(3) - v1aux(3) * v2aux(1))
    n_xyz(3) = v1aux(1) * v2aux(2) - v1aux(2) * v2aux(1)

    dd = n_xyz(1) * p_xyz(1,1) + n_xyz(2) * p_xyz(1,2) + n_xyz(3) * p_xyz(1,3)

    do c = 1, Grid % n_cells

      ! For every node:
      do n = 1, Grid % cells_n_nodes(c)

        res_dummy = n_xyz(1) * Grid % xn(Grid % cells_n(n,c))      &
                  + n_xyz(2) * Grid % yn(Grid % cells_n(n,c))      &
                  + n_xyz(3) * Grid % zn(Grid % cells_n(n,c))

        min_max_crit_1(c)= min(res_dummy, min_max_crit_1(c))
        min_max_crit_2(c)= max(res_dummy, min_max_crit_2(c))
      end do
    end do

    ! Simply interpolate linearly
    do c = 1, Grid % n_cells
      if (min_max_crit_1(c) < dd .and. min_max_crit_2(c) > dd) then
        prelim_vof(c) = 1.0 - (min_max_crit_2(c) - dd)  &
                      / (min_max_crit_2(c)-min_max_crit_1(c))
      else if (min_max_crit_2(c) <= dd) then
        prelim_vof(c) = 1.0
      end if
    end do

    ! Precision
    do c = 1, Grid % n_cells
      Vof % fun % n(c) = max(prelim_vof(c),Vof % fun % n(c))
    end do

  end do

  close(fu)

  call Work % Disconnect_Real_Cell(prelim_vof, min_max_crit_1, min_max_crit_2)

  end subroutine
