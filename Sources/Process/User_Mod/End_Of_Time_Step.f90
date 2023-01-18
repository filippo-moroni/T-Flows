!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
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
  
  double precision :: Cl = 0.0
  double precision :: Cd = 0.0
  
  double precision :: vis_dyn = 1.78e-5 ! Dynamic viscosity
   
  double precision :: n_versor(3)              ! Normal versor to the surface
  double precision :: t_versor(3)              ! Tangent versor to the surface
  double precision :: u_parallel 
  double precision :: tau_w
  
  double precision :: u_vel,v_vel,w_vel        ! Velocity components 
  double precision :: dsx, dsy, dsz, ds        ! Face element area vector components and total magnitude
  double precision :: wall_distance
   				
  integer :: c = 0
  
  logical :: a = .true.
  			
  character(len=1024) :: filename	
!==============================================================================!
  
  ! Consider only near-wall cells
  
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
  
  if (Grid % cell_near_wall(c) .eqv. a) then                        ! Near-wall cell selection
  
  	if (Grid % yc(c) .gt. -3 .and. Grid % yc(c) .lt. 3) then    ! Airfoil-wall selection
  	
  	u_vel = Flow % u % n(c)
  	v_vel = Flow % v % n(c)
  	w_vel = Flow % w % n(c) 
  	
  	ds =  Grid % s(c)
  	dsx = Grid % sx(c)
  	dsy = Grid % sy(c)
  	dsz = Grid % sz(c)
  	
  	wall_distance = Grid % wall_dist(c)
  	  	  	  	
  	n_versor(1) = dsx/ds                                      ! Components of the normal vector to the surface (wall face)                                           
  	n_versor(2) = dsy/ds
  	n_versor(3) = dsx/ds
  	
  	t_versor(1) =  n_versor(2)                                ! Tangent vector
  	t_versor(2) = -n_versor(1)
  	t_versor(3) =  n_versor(3)
  	
  	u_parallel = u_vel*t_versor(1) + v_vel*t_versor(2) + w_vel*t_versor(3)  ! Parallel velocity to the surface (first cell)
  	
  	tau_w = vis_dyn*u_parallel/wall_distance                    ! Wall shear stress
  	
  	Cl = Cl - Flow % p % n(c)*dsy - tau_w*dsy                   ! Lift
  	
  	Cd = Cd - Flow % p % n(c)*dsx + tau_w*dsx                   ! Drag
  	
  	end if
  
  end if
  	    
  end do
  
  ! Calculating the coefficients, being u_ref = 1 m/s, c = 1 m and b = 0.6 m:
  
  Cd = Cd / 0.36
  Cl = Cl / 0.36
  
  ! SubSnapshots creation
  
  if (this_proc < 10) write (filename, "(A2,I1,A4)") "Cl", this_proc, '.txt'			
  if (this_proc > 9)  write (filename, "(A2,I2,A4)") "Cl", this_proc, '.txt'
  
  if (this_proc < 10) write (filename, "(A2,I1,A4)") "Cd", this_proc, '.txt'			
  if (this_proc > 9)  write (filename, "(A2,I2,A4)") "Cd", this_proc, '.txt'
  
  ! Creation of the final output file with Cl and Cd 
  
  open(unit=911+this_proc,file = trim(filename),status='unknown',  &
       action='write',form='formatted',position="append")
  
    	write (911+this_proc, *) Cl, &
                             	 Cd						  
        	     	     
  close(911+this_proc)

  end subroutine
  
  
  
