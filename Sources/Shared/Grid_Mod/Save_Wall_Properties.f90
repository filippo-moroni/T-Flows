!==============================================================================!
  subroutine Save_Wall_Properties(Grid)
!------------------------------------------------------------------------------!
! Writes a .txt file that contains coordinates and wall faces properties:      !
! the normals to the surface, the face element area and the wall distance.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
!-----------------------------------[Locals]-----------------------------------!
    
  type Snapshot

  	real, allocatable    :: x_snap(:), y_snap(:), z_snap(:)     ! Coordinates
	real, allocatable    :: nx_snap(:), ny_snap(:), nz_snap(:)  ! Surface vector components
        real, allocatable    :: ds_snap(:)                          ! Surface vector magnitude
        real, allocatable    :: wd_snap(:)                          ! Wall distance

	real, allocatable    :: x_2(:), y_2(:), z_2(:)              ! Coordinates for the sorting

	integer, allocatable :: pos_ind(:)                          ! Position index

  end type

  type(Snapshot)  :: Snap
  
    real          :: nx,ny,nz  ! Normal vector components (versor).
    real          :: ds        ! Surface area of the element.     
    integer       :: s,c1,c2    
    integer       :: l = 0     ! Integer to count how many cells we have at the selected boundary.
    integer       :: c, iunit   
    character(SL) :: answer            
!==============================================================================!

  ! Display the BCs list and save the selected BC's properties 
  call Print_Bnd_Cond_List(Grid)
  print '(a)', ' #==================================='
  print '(a)', ' # Insert the BC name in upper cases '
  print '(a)', ' #-----------------------------------' 
  read *, answer
    
  ! Creation of a .txt file containing the cells' coordinates and their properties
  open (newunit=iunit, file='Wall_prop_uns.txt', form='formatted', status='unknown')
   
  do s = 1, Grid % n_faces
  
  	c1 = Grid % faces_c(1,s)
  	c2 = Grid % faces_c(2,s)
  	
  	! Take only faces at the boundaries  	
  	if(c2 < 0) then
  	
  		! Selection of the boundary we asked before 		
  		if(Grid % Bnd_Cond_Name(c2) .eq. answer) then
  		
  		ds = sqrt((Grid % sx(s))**2 + (Grid % sy(s))**2 + (Grid % sz(s))**2)
  		
  		nx = Grid % sx(s) / ds               ! Normal vector components
  		ny = Grid % sy(s) / ds
  		nz = Grid % sz(s) / ds
  	 	  	
		write (iunit, *) Grid % xc(c1),       &  ! x-coordinate
                     	         Grid % yc(c1),       &  ! y-coordinate
                    	         Grid % zc(c1),       &  ! z-coordinate                   
			         nx,                  &  ! x-component normal vector
			         ny,                  &  ! y-component normal vector
			         nz,                  &  ! z-component normal vector			     
                     	         ds,                  &  ! ds surface area
                    	         Grid % wall_dist(c1)    ! wall distance
                     
                                 l = l + 1               ! Parameter to count...
                                                         ! ... how many cells we have.     
  		end if
  	
  	end if
     
  end do
  
  close (iunit)

  ! Opening of the .txt file created to store the variables in the Snapshot

  allocate(Snap % x_snap(l));    Snap % x_snap = 0.0
  allocate(Snap % y_snap(l));    Snap % y_snap = 0.0
  allocate(Snap % z_snap(l));    Snap % z_snap = 0.0
  allocate(Snap % nx_snap(l));   Snap % nx_snap = 0.0
  allocate(Snap % ny_snap(l));   Snap % ny_snap = 0.0
  allocate(Snap % nz_snap(l));   Snap % nz_snap = 0.0
  allocate(Snap % ds_snap(l));   Snap % ds_snap = 0.0
  allocate(Snap % wd_snap(l));   Snap % wd_snap = 0.0
 
  open (newunit=iunit, file='Wall_prop_uns.txt', form='formatted', status='unknown')
   
  do c = 1, l
               
  	read (iunit,*)      Snap % x_snap(c),   &
                            Snap % y_snap(c),   &
                            Snap % z_snap(c),   &
                            Snap % nx_snap(c),  &
                            Snap % ny_snap(c),  &
                            Snap % nz_snap(c),  &
                            Snap % ds_snap(c),  &
                            Snap % wd_snap(c)

                                                                                                   			       
  end do
  
  close (iunit, status='delete')
  	
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

  open(newunit=iunit,file='Wall_Properties.txt',form='formatted',status='unknown')

  do c = 1, l

	write(iunit, *) Snap % x_snap(Snap % pos_ind(c)), ',', &
		        Snap % y_snap(Snap % pos_ind(c)), ',', &
		        Snap % z_snap(Snap % pos_ind(c)), ',', &
		        Snap % nx_snap(Snap % pos_ind(c)), ',', &
		        Snap % ny_snap(Snap % pos_ind(c)), ',', &
		        Snap % nz_snap(Snap % pos_ind(c)), ',', &
		        Snap % ds_snap(Snap % pos_ind(c)), ',', &
		        Snap % wd_snap(Snap % pos_ind(c))

  end do 
  
  close(iunit)

  deallocate(Snap % x_snap)     
  deallocate(Snap % y_snap)     
  deallocate(Snap % z_snap)    
  deallocate(Snap % nx_snap)   
  deallocate(Snap % ny_snap)   
  deallocate(Snap % nz_snap)   
  deallocate(Snap % ds_snap)   
  deallocate(Snap % wd_snap)   
  deallocate(Snap % pos_ind)
  deallocate(Snap % x_2)
  deallocate(Snap % y_2)
  deallocate(Snap % z_2)

  end subroutine
