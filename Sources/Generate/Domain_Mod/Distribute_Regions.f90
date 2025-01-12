!==============================================================================!
  subroutine Domain_Mod_Distribute_Regions(dom, Grid)
!------------------------------------------------------------------------------!
!   Distribute regions (defined in .dom file) to boundary conditions           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, i, j, k, n, c, r
  integer :: n_bnd                         ! number of boundary conditions
  integer :: is, js, ks, ie, je, ke, face
  integer :: ci, cj, ck
  logical :: found
!==============================================================================!

  !-----------------------------------------!
  !   Insertion of the boundary condition   !
  !        and materials information        !
  !-----------------------------------------!

  ! This is too much memory but that's OK 
  !  (+1 is to store the default values)
  allocate(Grid % bnd_cond % name(dom % n_regions + 1))

  ! Initialize number of boundary conditions
  n_bnd = 0

  do n = 1, dom % n_regions

    b = dom % regions(n) % block

    ! Block resolution
    ci = dom % blocks(b) % resolutions(1)-1
    cj = dom % blocks(b) % resolutions(2)-1
    ck = dom % blocks(b) % resolutions(3)-1

    ! Default values
    is = 1
    ie = ci
    js = 1
    je = cj
    ks = 1
    ke = ck

    ! Boundary conditions prescribed with mnemonics
    if(dom % regions(n) % face .eq. 'IMIN') then
      ie=1 
      face = 5
    else if(dom % regions(n) % face .eq. 'IMAX') then 
      is=ci
      face = 3
    else if(dom % regions(n) % face .eq. 'JMIN') then 
      je=1
      face = 2
    else if(dom % regions(n) % face .eq. 'JMAX') then 
      js=cj
      face = 4
    else if(dom % regions(n) % face .eq. 'KMIN') then 
      ke=1
      face = 1
    else if(dom % regions(n) % face .eq. 'KMAX') then 
      ks=ck
      face = 6

    ! Boundary conditions  prescribed explicitly
    !  (error prone and difficult, but might be usefull)
    else   
      is = dom % regions(n) % is
      js = dom % regions(n) % js
      ks = dom % regions(n) % ks
      ie = dom % regions(n) % ie
      je = dom % regions(n) % je
      ke = dom % regions(n) % ke
      face = 0
      if( (is .eq. ie).and.(is .eq.  1) ) face=5
      if( (is .eq. ie).and.(is .eq. ci) ) face=3
      if( (js .eq. je).and.(js .eq.  1) ) face=2
      if( (js .eq. je).and.(js .eq. cj) ) face=4
      if( (ks .eq. ke).and.(ks .eq.  1) ) face=1
      if( (ks .eq. ke).and.(ks .eq. ck) ) face=6
    end if

    ! Store boundary condition 
    if(face .ne. 0) then  

      found = .false. 
      do r=1,n_bnd
        if( Grid % bnd_cond % name(r) .eq.   &
            dom % regions(n) % name ) found = .true.
      end do
      if( .not. found) then
        n_bnd = n_bnd + 1
        Grid % bnd_cond % name(n_bnd) = dom % regions(n) % name
      end if

      do i=is,ie
        do j=js,je
          do k=ks,ke
            c = dom % blocks(b) % n_cells + (k-1)*ci*cj + (j-1)*ci + i   
            Grid % cells_c(face,c) = -n_bnd
          end do
        end do
      end do

    end if

  end do  !  n_regions

  ! Store the number of boundary conditions
  Grid % n_bnd_cond  = n_bnd

  call Grid % Print_Bnd_Cond_List()

  end subroutine
