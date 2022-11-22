!==============================================================================!
  subroutine Extrapolate_Normal_To_Front(Vof, Flow, phi, towards)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Field_Type),  target :: Flow
  type(Var_Type)            :: phi   ! it is t_0 or t_1
  integer, intent(in)       :: towards
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Front_Type), pointer :: Front
  real, pointer, contiguous :: phi_new(:), phi_old(:)
  real                      :: dx, dy, dz, nx, ny, nz, sx, sy, sz
  real                      :: dtau, flux_n, max_error
  integer                   :: t_iter, c, s, c1, c2, e1, e2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITER = 8
!==============================================================================!

  ! Take aliases
  Grid  => Flow % pnt_grid
  Front => Vof % Front

  !--------------------------------------------!
  !   Find dtau in coordinate invariant mode   !
  !--------------------------------------------!
  dtau = HUGE
  do s = 1, Grid % n_faces
    dx = abs(Grid % dx(s))
    dy = abs(Grid % dy(s))
    dz = abs(Grid % dz(s))
    dtau = min(dtau, sqrt((dx**2 + dy**2 + dz**2)))
  end do

  !---------------------------!
  !   Set initial condition   !
  !---------------------------!
  do c = 1, Grid % n_cells
    phi % x(c) = Flow % t % x(c)
    phi % y(c) = Flow % t % y(c)
    phi % z(c) = Flow % t % z(c)
  end do

  !--------------------------------!
  !   Pseudo time-step iteration   !
  !--------------------------------!
  do t_iter = 1, MAX_ITER

    ! The new becomes old
    do c = 1, Grid % n_cells
      phi_old(c) = phi_new(c)
    end do

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      ! old and bad: e1 = Front % elem_in_cell(c1)
      ! old and bad: e2 = Front % elem_in_cell(c2)
      e1 = Front % elems_at_face(1, s)
      e2 = Front % elems_at_face(2, s)

      if(e1 > 0 .or. e2 > 0) then

        ! Add the surfaces up (wouldn't bet my life this is a good idea)
        sx = 0.0
        sy = 0.0
        sz = 0.0
        if(e1 > 0) then
          sx = sx + Front % Elem(e1) % sx
          sy = sy + Front % Elem(e1) % sy
          sz = sz + Front % Elem(e1) % sz
        end if
        if(e2 > 0) then
          sx = sx + Front % Elem(e2) % sx
          sy = sy + Front % Elem(e2) % sy
          sz = sz + Front % Elem(e2) % sz
        end if

        ! Define element normals - this is the "new normal"
        nx = sx / sqrt(sx*sx + sy*sy + sz*sz)
        ny = sy / sqrt(sx*sx + sy*sy + sz*sz)
        nz = sz / sqrt(sx*sx + sy*sy + sz*sz)

        ! This is flux if the normal was in fact velocity.  By defintion (which
        ! we are trying to establish), it is oriented from vof=0 to vof=1)
        flux_n = nx * Grid % dx(s)  &
               + ny * Grid % dy(s)  &
               + nz * Grid % dz(s)

        ! If flux_n > 0, normal "flows" from c1 to c1, and also towards=1
        ! ensures we are extrapolating from vof=0 to vof=11
        if(flux_n > 0.0 .and. towards .eq. 1) then
          phi_new(c2) = phi_old(c2)                  &  ! unsteady term
                      - dtau * (  nx * phi % x(c1)   &  ! advection ...
                                + ny * phi % y(c1)   &  ! ... term ...
                                + nz * phi % z(c1))     ! ... is here
        else
          phi_new(c1) = phi_old(c1)                  &
                      + dtau * (  nx * phi % x(c2)   &
                                + ny * phi % y(c2)   &
                                + nz * phi % z(c2))

        end if

      end if  ! e1 > 0 .or. e2 > 0

    end do

    ! The new becomes old
    max_error = -HUGE
    do c = 1, Grid % n_cells
      max_error = max(max_error, abs(phi_old(c) - phi_new(c)))
    end do
    if(max_error < MICRO) goto 1

  end do  ! iter

1 continue

  end subroutine
