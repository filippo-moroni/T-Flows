!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
! This function is called at the end of time step and calculates Cd and Cl     !
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
    
  real :: vis_dyn = 1.78e-5        ! Dynamic viscosity
   
  real :: nx,ny                    ! Normal versor to the surface
  real :: tx,ty                    ! Tangent versor to the surface
  real :: u_parallel               ! Parallel velocity to the surface
  real :: tau_w                    ! Wall shear stress
  real :: pressure                 ! Pressure
  
  real :: u_vel,v_vel              ! Velocity components 
  real :: ds                       ! Area of the face at the boundary of an element
  real :: wall_distance            ! Wall distance          
  
  real :: Cl = 0.0
  real :: Cd = 0.0
  real :: Cl_tot = 0.0
  real :: Cd_tot = 0.0
  
  integer :: s,c1,c2 				
  integer :: j = 0
  integer :: time_interval = 20   ! This parameter controls how many ts are between the Cd and Cl saving
  integer :: b                    ! Integer to verify if it is time to save Cd and Cl
  
  logical :: exist
  			
  character(len=1024) :: filename	
!==============================================================================!

  ! Take aliases
 
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  
  ! Perform the drag and lift calculations if the time-step is the right one
  
  b = mod(n,time_interval)
  
  if (b .eq. 0) then
  
  ! Reset Cd and Cl values to zero
  
  Cd = 0.0
  Cl = 0.0 
  
  Cd_tot = 0.0
  Cl_tot = 0.0
    
  ! Browse through all cells and select the airfoil boundary
  
  do s = 1, Grid % n_faces
  
  	c1 = Grid % faces_c(1, s)
  	c2 = Grid % faces_c(2, s)
  	  	
  	! Take only faces at the boundaries

        if(c2 < 0) then
  	
  	! Then only pick faces at the boundary condition called AIRFOIL

        if(Grid % Bnd_Cond_Name(c2) .eq. 'AIRFOIL') then
  	
  		u_vel = Flow % u % n(c1)       ! Velocity components
  		v_vel = Flow % v % n(c1) 
  	
  		pressure = Flow % p % n(c1)    ! Pressure
  	
  		nx = Grid % sx(s)/Grid % s(s)  ! Normal vector to the surface
  		ny = Grid % sy(s)/Grid % s(s)
  		
  		tx =  ny                       ! Tangent vector to the surface
  		ty = -nx
  	
		ds = Grid % s(s)               ! Face element area 
  	
  		wall_distance = Grid % wall_dist(c1)         ! Wall distance
  	  	  	  	  	
  		u_parallel = u_vel*tx + v_vel*ty             ! Parallel velocity to the surface (first cell)
  	
  		tau_w = vis_dyn*u_parallel/wall_distance     ! Wall shear stress
  	
  		Cd = Cd - pressure*ds*nx + tau_w*ds*tx       ! Drag  
  		
  		Cl = Cl - pressure*ds*ny + tau_w*ds*ty       ! Lift
  	  	
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
  
