!==============================================================================!
  subroutine Vof_Interface_Cylinder(Vof,               &
                                    c,                 &
                                    p1_x, p1_y, p1_z,  &
                                    p2_x, p2_y, p2_z,  &
                                    radius, height,    &
                                    vof_int)
!------------------------------------------------------------------------------!
!   Computes volume fraction of cell at interface                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
  integer, intent(in)    :: c
  real                   :: p1_x, p1_y, p1_z
  real                   :: p2_x, p2_y, p2_z
  real                   :: radius, height
  real                   :: vof_int
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: N = 10000
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  logical                  :: l_cel
  integer                  :: nod, n_int, n_tot, n_tot_int, fu
  integer                  :: ee, n_cylinders, i_vari, j_vari, n_vari
  real                     :: r_num, res_func
  real                     :: xmin, xmax, ymin, ymax, zmin, zmax
  real                     :: p(3)
  real                     :: vof_0, vof_tol1, vof_tol2
  real                     :: avg_x, avg_y, avg_z
  real                     :: var_comb, var_comb_0, dist_cent
  real                     :: mean_x, mean_y, mean_z
  real                     :: points(N,3)
!==============================================================================!

  ! First take aliasesd
  Grid => Vof % pnt_grid

  ! Initialize variables
  var_comb   = HUGE
  var_comb_0 = 0.0
  vof_tol1   = MICRO
  vof_tol2   = MILI
  dist_cent  = Grid % vol(c) ** ONE_THIRD
  vof_int    = HUGE

  n_tot = 0
  n_int = 0

  ! find bounding box:
  xmin =  HUGE; ymin =  HUGE; zmin =  HUGE;
  xmax = -HUGE; ymax = -HUGE; zmax = -HUGE;

  do nod = 1, Grid % cells_n_nodes(c)
    xmin = min(xmin, Grid % xn(Grid % cells_n(nod,c)))
    ymin = min(ymin, Grid % yn(Grid % cells_n(nod,c)))
    zmin = min(zmin, Grid % zn(Grid % cells_n(nod,c)))
    xmax = max(xmax, Grid % xn(Grid % cells_n(nod,c)))
    ymax = max(ymax, Grid % yn(Grid % cells_n(nod,c)))
    zmax = max(zmax, Grid % zn(Grid % cells_n(nod,c)))
  end do

  l_cel = .false.

  do while ( (n_tot < N) .and. (abs(var_comb - var_comb_0) > vof_tol1 .or. &
             (var_comb / dist_cent) > vof_tol2) )

    l_cel = .false.

    ! check if p is inside cell
    do while (l_cel .eqv. .false.)
      call random_number(r_num)
      p(1) = xmin + (xmax-xmin) * r_num
      call random_number(r_num)
      p(2) = ymin+ (ymax-ymin) * r_num
      call random_number(r_num)
      p(3) = zmin+ (zmax-zmin) * r_num
      l_cel = Grid % Is_Point_In_Cell(c, p(1), p(2), p(3))
    end do
    n_tot = n_tot + 1

    points(n_tot,:) = p(:)

    ! Check if p is inside function:
    res_func = ( ( (p2_y-p1_y)*(p1_z-p(3))        &
                  -(p1_y-p(2))*(p2_z-p1_z))**2    &
                +( (p2_x-p1_x)*(p1_z-p(3))        &
                  -(p1_x-p(1))*(p2_z-p1_z))**2    &
                +( (p2_x-p1_x)*(p1_y-p(2))        &
                  -(p1_x-p(1))*(p2_y-p1_y))**2 )  &
                / (radius*height) ** 2
    if (res_func <= 1.0) then
      n_int = n_int + 1
    end if

    ! Compute maximum difference c/r center:
    if (mod(n_tot,500) == 0 .and. n_tot > 0) then
      var_comb_0 = var_comb
      vof_int = real(n_int) / real(n_tot)
      mean_x = sum(points(1: n_tot,1)) / real(n_tot)
      mean_y = sum(points(1: n_tot,2)) / real(n_tot)
      mean_z = sum(points(1: n_tot,3)) / real(n_tot)
      var_comb = abs(mean_x-Grid % xc(c))   &
               + abs(mean_y-Grid % yc(c))   &
               + abs(mean_z-Grid % zc(c))
    end if
  end do

  end subroutine
