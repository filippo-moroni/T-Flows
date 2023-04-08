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
!==============================================================================!

  integer           :: tot_iter = 20000
  integer           :: l
  real, allocatable :: Cl(:) = 0.0
  real, allocatable :: Cd(:) = 0.0
  
  real, allocatable :: Cl_tot(:) = 0.0
  real, allocatable :: Cd_tot(:) = 0.0
  
  logical :: exist		
  character(len=1024) :: filename	
  
! Creation of the total output of Cd and Cl
  if (this_proc .eq. n_proc) then
  
  allocate

  do j = 1, n_proc									             
    
      if (j < 10)               write (filename, "(A5,I1,A4)") "Cd+Cl", j, '.txt'				     
      if (j > 9 .and. j < 100)  write (filename, "(A5,I2,A4)") "Cd+Cl", j, '.txt'
      if (j > 99)               write (filename, "(A5,I3,A4)") "Cd+Cl", j, '.txt'
  
      open(unit=651+j,file = trim(filename),form='formatted',status='unknown')	
      
      read(651+j, *, end=10) Cd, Cl
	
      10 close(unit=651+j,status='delete')
      
      Cd_tot = Cd_tot + Cd_read 	            
      Cl_tot = Cl_tot + Cl_read  
           
  end do
  
  inquire(file="Cd+Cl_tot.txt", exist=exist)
  
  if (exist) then
    open(unit=796, file="Cd+Cl_tot.txt", status='old', position='append', action='write', form='formatted')
  else
    open(unit=796, file="Cd+Cl_tot.txt", status='new', action='write', form='formatted')
  end if
         
  write(796, *) n, &
  	        Cd_tot, &
                Cl_tot
  		
  close(796)
  
  end if

  end subroutine
