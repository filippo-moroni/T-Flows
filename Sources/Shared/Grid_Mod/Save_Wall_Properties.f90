!==============================================================================!
  subroutine Save_Wall_Properties(Grid)
!------------------------------------------------------------------------------!
! Writes a .txt file that contains coordinates and wall faces properties:      !
! the 3 surface vector components, the face area and the wall distance         !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
!-----------------------------------[Locals]-----------------------------------!
    
  type Snapshot

  	real, allocatable :: x_snap(:), y_snap(:), z_snap(:)            ! Coordinates
	real, allocatable :: dsx_snap(:), dsy_snap(:), dsz_snap(:)      ! Surface vector components
        real, allocatable :: ds_snap(:)                                 ! Surface vector magnitude
        real, allocatable :: wd_snap(:)                                 ! Wall distance

	real, allocatable :: x_2(:), y_2(:), z_2(:)

	integer, allocatable :: pos_ind(:)

  end type

  type(Snapshot) :: Snap
      
    integer       :: c = 0 
    integer       :: l = 0
    integer       :: a = 0
    integer       :: answer
    
    integer       :: n_bnd_sect
     
    !character       :: answer
             
!==============================================================================!
 
  print '(a)', ' #=================================='
  print '(a)', ' # Insert the BC n° '
  print '(a)', ' #----------------------------------'
  
  call Print_Bnd_Cond_List(Grid)
  
  read *, answer
    
  ! Creation of a .txt file containing the cells' coordinates and their properties

  open (unit=1, file='Wall_prop_uns.txt', form='formatted', status='unknown')
   
  do c = 1, Grid % n_cells
  
  do a = 1, Grid % n_bnd_cond
  
  if (Grid % cells_bnd_color(c,a) .eq. answer) then       ! how to open only the airfoil boundary?
  
  !if (Grid % bnd_cond % name(a) == answer) then
  	
	write (1, *) Grid % xc(c),       & ! x-coordinate
                     Grid % yc(c),       & ! y-coordinate
                     Grid % zc(c),       & ! z-coordinate                   
		     Grid % sx(c),       & ! dsx surface vector component
  		     Grid % sy(c),       & ! dsy surface vector component
	             Grid % sz(c),       & ! dsz surface vector component
                     Grid % s(c),        & ! ds surface vector magnitude
                     Grid % wall_dist(c)   ! wall distance
                     
                     l = l + 1
        
  end if
  
  end do
      
  end do
  
  close (1)

  ! Opening of the .txt file created to store the variables in the Snapshot

  allocate(Snap % x_snap(l));     Snap % x_snap = 0.0
  allocate(Snap % y_snap(l));     Snap % y_snap = 0.0
  allocate(Snap % z_snap(l));     Snap % z_snap = 0.0
  allocate(Snap % dsx_snap(l));   Snap % dsx_snap = 0.0
  allocate(Snap % dsy_snap(l));   Snap % dsy_snap = 0.0
  allocate(Snap % dsz_snap(l));   Snap % dsz_snap = 0.0
  allocate(Snap % ds_snap(l));    Snap % ds_snap = 0.0
  allocate(Snap % wd_snap(l));    Snap % wd_snap = 0.0
 
  open (unit=2, file='Wall_prop_uns.txt', form='formatted', status='unknown')
   
  do c = 1, l
               
  	read (2,*)          Snap % x_snap(c),   &
                            Snap % y_snap(c),   &
                            Snap % z_snap(c),   &
                            Snap % dsx_snap(c), &
                            Snap % dsy_snap(c), &
                            Snap % dsz_snap(c), &
                            Snap % ds_snap(c),  &
                            Snap % wd_snap(c)

                                                                                                   			       
  end do
  
  close (2, status='delete')
  	
  ! Ordering the cells through their coordinates

  allocate(Snap % pos_ind(l));  Snap % pos_ind = 0
  allocate(Snap % x_2(l));      Snap % x_2 = 0.0
  allocate(Snap % y_2(l));      Snap % y_2 = 0.0
  allocate(Snap % z_2(l));      Snap % z_2 = 0.0

  Snap % x_2(:) = HUGE
  Snap % y_2(:) = HUGE

  do c = 1, l 

  	Snap % pos_ind(c) = c
        Snap % x_2(c) = min(Snap % x_2(c), Snap % x_snap(c))
        Snap % y_2(c) = min(Snap % y_2(c), Snap % y_snap(c))
	Snap % z_2(c) = Snap % z_snap(c)

  end do

  Call Sort % Three_Real_Carry_Int(Snap % x_2, Snap % y_2, Snap % z_2, Snap % pos_ind)

  ! Writing the file with ordered cells according to their coordinates

  open(unit=3,file='Wall_Properties.txt',form='formatted',status='unknown')

  do c = 1, l

	write(3, *) Snap % x_snap(Snap % pos_ind(c)), ',', &
		    Snap % y_snap(Snap % pos_ind(c)), ',', &
		    Snap % z_snap(Snap % pos_ind(c)), ',', &
		    Snap % dsx_snap(Snap % pos_ind(c)), ',', &
		    Snap % dsy_snap(Snap % pos_ind(c)), ',', &
		    Snap % dsz_snap(Snap % pos_ind(c)), ',', &
		    Snap % ds_snap(Snap % pos_ind(c)), ',', &
		    Snap % wd_snap(Snap % pos_ind(c))

  end do 
  
  close(3)

  deallocate(Snap % x_snap)     
  deallocate(Snap % y_snap)     
  deallocate(Snap % z_snap)    
  deallocate(Snap % dsx_snap)   
  deallocate(Snap % dsy_snap)   
  deallocate(Snap % dsz_snap)   
  deallocate(Snap % ds_snap)   
  deallocate(Snap % wd_snap)   
  deallocate(Snap % pos_ind)
  deallocate(Snap % x_2)
  deallocate(Snap % y_2)
  deallocate(Snap % z_2)

  end subroutine






