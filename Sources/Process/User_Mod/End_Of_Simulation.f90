!==============================================================================!
  subroutine User_Mod_End_Of_Simulation(Flow, Turb, Vof, Swarm, curr_dt, time)
!------------------------------------------------------------------------------!
! This function is called at the end of simulation and recollects              !
! the Cd and Cl SubSnapshots savings.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  integer              :: tot_iter = 15000		! Maximum number of iterations in a single run, to allocate memory.
  integer              :: c = 0
  integer              :: i, j
  real,    allocatable :: Cl(:) = 0.0
  real,    allocatable :: Cd(:) = 0.0
  integer, allocatable :: n(:) = 0
  
  real,    allocatable :: Cl_tot(:) = 0.0
  real,    allocatable :: Cd_tot(:) = 0.0
  integer, allocatable :: n_tot(:) = 0
   
  logical :: exist		
  character(len=1024) :: filename
!==============================================================================!
	  
! Creation of the total output of Cd and Cl
  if (this_proc .eq. n_proc) then
  
  allocate(n(tot_iter));	n(:)=0
  allocate(Cd(tot_iter));	Cd(:)=0.0
  allocate(Cl(tot_iter));	Cl(:)=0.0
  
  allocate(n_tot(tot_iter));	n_tot(:)=0
  allocate(Cd_tot(tot_iter));	Cd_tot(:)=0.0
  allocate(Cl_tot(tot_iter));	Cl_tot(:)=0.0
  
  ! Open and sum all the contributions from the SubSnapshot files
  do j = 1, n_proc
  
      c = 0
    
      if (j < 10)               write (filename, "(A5,I1,A4)") "Cd+Cl", j, '.txt'				     
      if (j > 9 .and. j < 100)  write (filename, "(A5,I2,A4)") "Cd+Cl", j, '.txt'
      if (j > 99)               write (filename, "(A5,I3,A4)") "Cd+Cl", j, '.txt'
  
      open(unit=651+j,file = trim(filename),form='formatted',status='old')	
      
      do
      
      read(651+j, *, end=17) n(c+1), Cd(c+1), Cl(c+1)
      
      n_tot(c+1) = n_tot(c+1) + n(c+1)
      Cd_tot(c+1) = Cd_tot(c+1) + Cd(c+1)
      Cl_tot(c+1) = Cl_tot(c+1) + Cl(c+1)
   
      c = c + 1		! Total number of rows in the SubSnapshot files
      
      end do
	
      17 close(unit=651+j,status='delete')  
    
  end do
  
  ! Debugging through the check of the number of iterations
  n_tot = n_tot/c
  
  ! Create the final overall file or update it
  inquire(file="Cd+Cl_tot.txt", exist=exist)
  
  if (exist) then
    open(unit=796, file="Cd+Cl_tot.txt", status='old', position='append', action='write', form='formatted')
  else
    open(unit=796, file="Cd+Cl_tot.txt", status='new', action='write', form='formatted')
  end if
  
  do i = 1, c
 	 write(796, *) n_tot(i), &
  	               Cd_tot(i), &
               	       Cl_tot(i)	
 	 close(796)
	 
  end do
  
  deallocate(n)
  deallocate(Cd)
  deallocate(Cl)
  deallocate(n_tot)
  deallocate(Cd_tot)
  deallocate(Cl_tot)
  
  end if ! Main if to run with just one processor

  end subroutine
