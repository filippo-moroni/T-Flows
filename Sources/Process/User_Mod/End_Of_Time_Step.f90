!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
! This function is called at the end of time step and calculates Cd and Cl.    !
! every 20 time steps.                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: n         ! time step
  integer, intent(in)         :: n_stat_t  ! start time step for Turb. stat.
  integer, intent(in)         :: n_stat_p  ! start time step for Swarm. stat.
  real,    intent(in)         :: time      ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
    
  real :: vis_dyn = 1.78e-5 ! Dynamic viscosity
   
  real :: n_versor(3)       ! Normal versor to the surface
  real :: t_versor(3)       ! Tangent versor to the surface
  real :: u_parallel 
  real :: tau_w
  real :: pressure
  
  real :: u_vel,v_vel,w_vel       ! Velocity components 
  real :: dsx, dsy, dsz, ds       ! Face element area vector components and total magnitude
  real :: wall_distance
  
  real :: Cl = 0.0
  real :: Cd = 0.0
  real :: Cl_tot = 0.0
  real :: Cd_tot = 0.0
   				
  integer :: c = 0
  integer :: j = 0
  integer :: time_interval = 20   ! This parameter controls how many ts are between the Cd and Cl saving
  integer :: b
  
  logical :: a = .true.
  logical :: exist
  			
  character(len=1024) :: filename	
!==============================================================================!

  ! Take aliases
 
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  
  ! Perform the drag and lift calculation if the ts is the right one
  
  b = mod(n,time_interval)
  
  if (b .eq. 0) then
  
  ! Reset Cd and Cl values to zero
  
  Cd = 0.0
  Cl = 0.0 
  
  Cd_tot = 0.0
  Cl_tot = 0.0
    
  ! Consider only near-wall cells
  
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
  
  if (Grid % cell_near_wall(c) .eqv. a) then                           ! Near-wall cell selection
  
  	if (Grid % yc(c) .gt. -2.0 .and. Grid % yc(c) .lt. 2.0) then   ! Airfoil-wall selection
  	
  	u_vel = Flow % u % n(c)
  	v_vel = Flow % v % n(c)
  	w_vel = Flow % w % n(c) 
  	
  	pressure = Flow % p % n(c)
  	
  	dsx = Grid % sx(c)
  	dsy = Grid % sy(c)
  	dsz = Grid % sz(c)
  	
  	ds = sqrt(dsx**2 + dsy**2 + dsz**2)
  	
  	wall_distance = Grid % wall_dist(c)
  	  	  	  	
  	n_versor(1) = dsx/ds                 ! Components of the normal vector to the surface (wall face)                                           
  	n_versor(2) = dsy/ds
  	n_versor(3) = dsz/ds
  	
  	t_versor(1) =  n_versor(2)           ! Tangent vector
  	t_versor(2) = -n_versor(1)
  	t_versor(3) =  n_versor(3)
  	
  	u_parallel = u_vel*t_versor(1) + v_vel*t_versor(2) + w_vel*t_versor(3)  ! Parallel velocity to the surface (first cell)
  	
  	tau_w = vis_dyn*u_parallel/wall_distance                  ! Wall shear stress
  	
  	Cd = Cd - pressure*ds*n_versor(1) + tau_w*ds*t_versor(1)  ! Drag  
  		
  	Cl = Cl - pressure*ds*n_versor(2) + tau_w*ds*t_versor(2)  ! Lift
  	  	
  	end if
  
  end if
  	    
  end do
  
  ! Calculating the coefficients and inverting the signs
  ! Reference quantities: u_ref = 1 m/s, c = 1 m and b = 0.6 m
  
  Cd = -Cd / 0.36
  
  Cl = -Cl / 0.36
  
  ! SubSnapshots creation
  
  if (this_proc < 10) write (filename, "(A5,I1,A4)") "Cd+Cl", this_proc, '.txt'			
  if (this_proc > 9)  write (filename, "(A5,I2,A4)") "Cd+Cl", this_proc, '.txt'
    
  ! Creation of the SubSnapshots with Cd and Cl 
  
  open(unit=911+this_proc,file = trim(filename),status='unknown',form='formatted')
  
    	write (911+this_proc, *) Cd, &
                             	 Cl						       	     	     
  close(911+this_proc)
  
!----------------------------------------------------------!

  ! Creation of the total output
  
  if (this_proc .eq. n_proc) then
  
  do j = 1, n_proc									             
    
      if (j < 10) write (filename, "(A5,I1,A4)") "Cd+Cl", j, '.txt'				     
      if (j > 9)  write (filename, "(A5,I2,A4)") "Cd+Cl", j, '.txt'				     
  
      open(unit=606+j,file = trim(filename),form='formatted',status='unknown')	
      
      read (606+j, *,end=10)  Cd, &				             
		              Cl
		              
      10 close(unit=606+j, status='delete')
	
      Cd_tot = Cd_tot + Cd
      	            
      Cl_tot = Cl_tot + Cl
           
  end do
  
  inquire(file="Cd+Cl_tot.txt", exist=exist)
  
  if (exist) then
    open(201, file="Cd+Cl_tot.txt", status="old", position="append", action="write")
  else
    open(201, file="Cd+Cl_tot.txt", status="new", action="write")
  end if
         
  write(201, *) n, &
  		Cd_tot, &
  		Cl_tot
  		
  close(201)
  
  end if
  
  end if ! Main if
   
  end subroutine
  
  
  
