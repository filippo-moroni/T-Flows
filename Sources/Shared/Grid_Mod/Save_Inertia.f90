!==============================================================================!
  subroutine Save_Inertia(Grid)
!------------------------------------------------------------------------------!
! Writes a .txt file that contains the inertia moments and the cells volume,   !
! used for post-processing of the TVM model.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
!-----------------------------------[Locals]-----------------------------------!
    
  type Snapshot

  	real, allocatable    :: x_snap(:), y_snap(:), z_snap(:)        ! Coordinates
	real, allocatable    :: ixx_snap(:), iyy_snap(:), izz_snap(:)  ! Inertia moments
        real, allocatable    :: ixy_snap(:), ixz_snap(:), iyz_snap(:)
        real, allocatable    :: cv_snap(:)                             ! Cell volume
	real, allocatable    :: x_2(:), y_2(:), z_2(:)                 ! Coordinates for sorting                    
	integer, allocatable :: pos_ind(:)                             ! Integer for sorting

  end type

  type(Snapshot) :: Snap
      
    integer :: c = 0 
    integer :: iunit          
!==============================================================================!

  ! Creation of a .txt file containing the cells' coordinates and the inertia moments
  open (newunit=iunit, file='Inertia_uns.txt', form='formatted', status='unknown')

  do c = 1, Grid % n_cells 
  
	write (iunit, *) Grid % xc(c), &
                         Grid % yc(c), &
                         Grid % zc(c), &                    
		         Grid % ixx(c), &
  		         Grid % iyy(c), &
	                 Grid % izz(c), &
                         Grid % ixy(c), &
		         Grid % ixz(c), &
                         Grid % iyz(c), &
                         Grid % vol(c)                                         
  end do
  
  close (iunit)

  ! Opening of the .txt file created to store the variables in the Snapshot
  allocate(Snap % x_snap(Grid % n_cells));     Snap % x_snap = 0.0
  allocate(Snap % y_snap(Grid % n_cells));     Snap % y_snap = 0.0
  allocate(Snap % z_snap(Grid % n_cells));     Snap % z_snap = 0.0
  allocate(Snap % ixx_snap(Grid % n_cells));   Snap % ixx_snap = 0.0
  allocate(Snap % iyy_snap(Grid % n_cells));   Snap % iyy_snap = 0.0
  allocate(Snap % izz_snap(Grid % n_cells));   Snap % izz_snap = 0.0
  allocate(Snap % ixy_snap(Grid % n_cells));   Snap % ixy_snap = 0.0
  allocate(Snap % ixz_snap(Grid % n_cells));   Snap % ixz_snap = 0.0
  allocate(Snap % iyz_snap(Grid % n_cells));   Snap % iyz_snap = 0.0
  allocate(Snap % cv_snap(Grid % n_cells));    Snap % cv_snap = 0.0
  
  open (newunit=iunit, file='Inertia_uns.txt', form='formatted', status='unknown')

  do c = 1, Grid % n_cells
               
  	read (iunit,*)      Snap % x_snap(c), &
                            Snap % y_snap(c), &
                            Snap % z_snap(c), &
                            Snap % ixx_snap(c), &
                            Snap % iyy_snap(c), &
                            Snap % izz_snap(c), &
                            Snap % ixy_snap(c), &
                            Snap % ixz_snap(c), &
                            Snap % iyz_snap(c), &
                            Snap % cv_snap(c)                                                                                                  			       
  end do
  
  close (iunit, status='delete')
  	
  ! Ordering the cells through their coordinates
  allocate(Snap % pos_ind(Grid % n_cells));  Snap % pos_ind = 0
  allocate(Snap % x_2(Grid % n_cells));      Snap % x_2 = 0.0
  allocate(Snap % y_2(Grid % n_cells));      Snap % y_2 = 0.0
  allocate(Snap % z_2(Grid % n_cells));      Snap % z_2 = 0.0

  Snap % x_2(:) = HUGE
  Snap % y_2(:) = HUGE

  do c = 1, Grid % n_cells 

  	Snap % pos_ind(c) = c
        Snap % x_2(c) = min(Snap % x_2(c), Snap % x_snap(c))
        Snap % y_2(c) = min(Snap % y_2(c), Snap % y_snap(c))
	Snap % z_2(c) = Snap % z_snap(c)

  end do

  Call Sort % Three_Real_Carry_Int(Snap % x_2, Snap % y_2, Snap % z_2, Snap % pos_ind)

  ! Writing the file with ordered cells according to their coordinates
  open(newunit=iunit,file='Inertia.txt',form='formatted',status='unknown')

  do c = 1, Grid % n_cells

	write(iunit, *) Snap % x_snap(Snap % pos_ind(c)), ',', &
		        Snap % y_snap(Snap % pos_ind(c)), ',', &
		        Snap % z_snap(Snap % pos_ind(c)), ',', &
		        Snap % ixx_snap(Snap % pos_ind(c)), ',', &
		        Snap % iyy_snap(Snap % pos_ind(c)), ',', &
		        Snap % izz_snap(Snap % pos_ind(c)), ',', &
		        Snap % ixy_snap(Snap % pos_ind(c)), ',', &
		        Snap % ixz_snap(Snap % pos_ind(c)), ',', &
		        Snap % iyz_snap(Snap % pos_ind(c)), ',', &
		        Snap % cv_snap(Snap % pos_ind(c))
  end do 
  
  close(iunit)

  deallocate(Snap % x_snap)     
  deallocate(Snap % y_snap)     
  deallocate(Snap % z_snap)    
  deallocate(Snap % ixx_snap)   
  deallocate(Snap % iyy_snap)   
  deallocate(Snap % izz_snap)   
  deallocate(Snap % ixy_snap)   
  deallocate(Snap % ixz_snap)   
  deallocate(Snap % iyz_snap)
  deallocate(Snap % cv_snap)
  deallocate(Snap % pos_ind)
  deallocate(Snap % x_2)
  deallocate(Snap % y_2)
  deallocate(Snap % z_2)

  end subroutine
