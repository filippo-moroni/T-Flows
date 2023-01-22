!==============================================================================!
  subroutine User_Mod_End_Of_Simulation(Flow, Turb, Vof, Swarm, &
                                        curr_dt, time, domain)
!------------------------------------------------------------------------------!
! This function is called at the end of the simulation and uses Process        !
! to calculate the wall distance and storing wall cells properties (dsi).      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
  integer,           optional :: domain
!==============================================================================!

  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w

  type Snapshot

  real, allocatable :: x_snap(:), y_snap(:), z_snap(:)            ! Coordinates
	real, allocatable :: dsx_snap(:), dsy_snap(:), dsz_snap(:)      ! Surface vector components
        real, allocatable :: ds_snap(:)                           ! Surface vector magnitude
        real, allocatable :: wd_snap(:)                           ! Wall distance

	real, allocatable :: x_2(:), y_2(:), z_2(:)

	integer, allocatable :: pos_ind(:)

  end type

  type(Snapshot) :: Snap
      
    integer       :: c = 0 
    integer       :: j = 0 
    integer       :: n_of_cells = 10000000 
    integer       :: n_c_tot = 0
      
    logical       :: a = .true.
    logical       :: exist
    
    real          :: ds_wall
    
    character(len=1024) :: filename	
                 
!==============================================================================!

  ! Take aliases
 
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Check if the 'Wall_Properties.txt' file already exists
  
  inquire(file="Wall_Properties.txt", exist=exist)
  
  if (exist) then
    
  else
  
  ! SubSnapshots creation
  
  if (this_proc < 10) write (filename, "(A9,I1,A4)") "Wall_Prop", this_proc, '.txt'			
  if (this_proc > 9)  write (filename, "(A9,I2,A4)") "Wall_Prop", this_proc, '.txt'
    
  ! Creation of the SubSnapshots with the wall properties and coordinates
  
  open(unit=911+this_proc,file = trim(filename),status='unknown',form='formatted')
  
  do c = 1, Grid % n_cells 
  
  if (Grid % cell_near_wall(c) .eqv. a) then                       ! Near-wall cell selection
  
  	if (Grid % yc(c) .gt. -2.0 .and. Grid % yc(c) .lt. 2.0) then   ! Airfoil-wall selection
  	
  	ds_wall = sqrt((Grid % sx(c))**2 + (Grid % sy(c))**2 + (Grid % sz(c))**2)
  	
	  write (911+this_proc, *) Grid % xc(c),       & ! x-coordinate
                             Grid % yc(c),       & ! y-coordinate
                     		     Grid % zc(c),       & ! z-coordinate                   
		    		                 Grid % sx(c),       & ! dsx surface vector component
  		    		               Grid % sy(c),       & ! dsy surface vector component
	             		           Grid % sz(c),       & ! dsz surface vector component
	             		           ds_wall,            & ! ds surface vector magnitude
                     		     Grid % wall_dist(c)   ! wall distance
                    
        end if
        
  end if
  
  end do 	
                             							       	     	     
  close(911+this_proc)
  
  ! Creation of the total output
  
  if (this_proc .eq. n_proc) then
 
  allocate(Snap % x_snap(n_of_cells));     Snap % x_snap = 0.0
  allocate(Snap % y_snap(n_of_cells));     Snap % y_snap = 0.0
  allocate(Snap % z_snap(n_of_cells));     Snap % z_snap = 0.0
  allocate(Snap % dsx_snap(n_of_cells));   Snap % dsx_snap = 0.0
  allocate(Snap % dsy_snap(n_of_cells));   Snap % dsy_snap = 0.0
  allocate(Snap % dsz_snap(n_of_cells));   Snap % dsz_snap = 0.0
  allocate(Snap % ds_snap(n_of_cells));    Snap % ds_snap = 0.0
  allocate(Snap % wd_snap(n_of_cells));    Snap % wd_snap = 0.0
  
  do j = 1, n_proc									             
    
      if (j < 10) write (filename, "(A9,I1,A4)") "Wall_Prop", j, '.txt'			
      if (j > 9)  write (filename, "(A9,I2,A4)") "Wall_Prop", j, '.txt'				     
  
      open(unit=606+j,file = trim(filename),form='formatted',status='unknown')	
      
      do
      
      read (606+j, *,end=10)  Snap % x_snap(n_c_tot+1), &
                              Snap % y_snap(n_c_tot+1), &
                              Snap % z_snap(n_c_tot+1), &
                              Snap % dsx_snap(n_c_tot+1), &
                              Snap % dsy_snap(n_c_tot+1), &
                              Snap % dsz_snap(n_c_tot+1), &
                              Snap % ds_snap(n_c_tot+1), &
                              Snap % wd_snap(n_c_tot+1)				             
		              
		              n_c_tot = n_c_tot + 1
      end do
		              		              
      10 close(unit=606+j, status='delete')
           
  end do
  	
  ! Ordering the cells through their coordinates

  allocate(Snap % pos_ind(n_c_tot));  Snap % pos_ind = 0
  allocate(Snap % x_2(n_c_tot));      Snap % x_2 = 0.0
  allocate(Snap % y_2(n_c_tot));      Snap % y_2 = 0.0
  allocate(Snap % z_2(n_c_tot));      Snap % z_2 = 0.0

  Snap % x_2(:) = HUGE
  Snap % y_2(:) = HUGE

  do c = 1, n_c_tot 

  	Snap % pos_ind(c) = c
        Snap % x_2(c) = min(Snap % x_2(c), Snap % x_snap(c))
        Snap % y_2(c) = min(Snap % y_2(c), Snap % y_snap(c))
	Snap % z_2(c) = Snap % z_snap(c)

  end do

  Call Sort % Three_Real_Carry_Int(Snap % x_2, Snap % y_2, Snap % z_2, Snap % pos_ind)

  ! Writing the file with ordered cells according to their coordinates

  open(unit=3,file='Wall_Properties.txt',form='formatted',status='unknown')

  do c = 1, n_c_tot

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
  
  n_c_tot = 0
  
  end if
  
  end if ! Main if of inquire
  
  end subroutine
