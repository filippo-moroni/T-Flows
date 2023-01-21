!==============================================================================!
  program Calculate_Total_Average
!------------------------------------------------------------------------------!
!          This program calculates the total-average of the                    !       
!                            Snapshots                                         !                                                                
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
      
  ! Variables definition
  
  character(len=1024) :: filename
  
  integer :: n_of_cells = 10000000                     ! Rough estimation of the number of cells of the grid to initialize the matrices. 
                                                       ! It must be greater than the real value of cells.
                                                       
  integer :: n_c_tot = 0                               ! Real total number of cells calculated during the 'do' cycle.
  integer :: n_c_tot_red = 0                           ! Reduced number of cells once we perform the spanwise average.                               
  
  integer :: i = 0
  integer :: j = 0
  
  integer :: n_layers                                  ! Number of layers in z-direction. This parameter is asked to the user.
  integer :: a = 1                                     ! Parameter to make the 'do' cycle work.
   
  integer :: init_ts                                   ! Corresponding time-step at which we start to calculate the time average.
  integer :: end_ts                                    ! Corresponding time-step at which we stop to calculate the time average.
  integer :: ts_interval                               ! Interval between each time step.
  integer :: n_snaps                                   ! Number of Snapshots.
  integer :: answer                                    ! Integer to consider the turbulent viscosity or not.
  integer :: ts                                        ! Extra variable to perform the Reynolds stress calculations.
  integer :: selection 		                       ! Integer to choose between performing the total average of some Snapshots 
                                                       ! or only the spanwise average of a single Snapshot.
    
  ! Definition of some matrices and variables
 
  double precision, allocatable :: snap(:,:)           ! To store the Snapshots          
  double precision, allocatable :: sum_matrix(:,:)     ! To store the summations in time
  double precision, allocatable :: sum_matrix_2(:,:)   ! To store the summations in spanwise direction
  
  double precision, allocatable :: inertia(:,:)        ! To store inertia (TVM model)
  double precision :: tau(3,3)                         ! To store subgrid stresses for TVM model at each iteration
  double precision :: ip(6)                            ! To store the inertia prime matrix for the TVM model at each iteration
  double precision :: sgv(3,3)                         ! To store the subgrid viscosity tensor
  double precision :: cv                               ! To store the cell volume
  double precision :: sij(3,3)                         ! To store the rate of deformation tensor
  double precision :: ni_sgs                           ! To store the sub-grid viscosity for Dynamic Smagorinsky simulations
  
  double precision, allocatable :: TKE_vector(:)       ! Turbulent kinetic energy column vector
  
  ! Definition of the matrices and vector that store the fluctuations
  
  double precision, allocatable :: u_prime(:,:)         ! u' fluctuations
  double precision, allocatable :: v_prime(:,:)         ! v' fluctuations
  double precision, allocatable :: w_prime(:,:)         ! w' fluctuations
  
  double precision, allocatable :: uu(:,:)              ! Instantaneous Reynolds stress u'u'
  double precision, allocatable :: vv(:,:)              ! Instantaneous Reynolds stress v'v'
  double precision, allocatable :: ww(:,:)              ! Instantaneous Reynolds stress w'w'
  double precision, allocatable :: rs(:,:)              ! Instantaneous Reynolds stress u'v'  
  
  double precision, allocatable :: sum_uu(:)            ! Sum of instantaneous Reynolds stress u'u'
  double precision, allocatable :: sum_vv(:)            ! Sum of instantaneous Reynolds stress v'v'
  double precision, allocatable :: sum_ww(:)            ! Sum of instantaneous Reynolds stress w'w'
  double precision, allocatable :: sum_rs(:)            ! Sum of instantaneous Reynolds stress u'v'
    
!------------------------------------------------------------------------------!
  ! Display the logo
   
  call logo()
!------------------------------------------------------------------------------!

  ! Ask the user what kind of file we are averaging
  
  print *,'What kind of average do you want to perform?'
  print *,'Press 1 to perform the time and spanwise averages of some Snapshots.'
  print *,'Press 2 to perform only the spanwise average of a specific Snapshot.'
  read *, selection
  
  ! Ask the user to insert the number of layers
  
  print *,'Insert the number of layers in z-direction.'
  read *, n_layers
  
  ! LES model used
  
  print *,'What LES model did you use?'
  print *,'Press 1 for Dynamic Smagorinsky (DSM)'
  print *,'Press 2 for Tensorial-Viscosity model (TVM).' 
  print *,'Press 3 for Implicit LES (ILES).' 
  read *, answer 
  
  !-----------------------------------------------------------------------------!
  
  if (selection .eq. 1) then ! Total average treatment
  
  ! Ask the user to insert the initial time-step, the final time-step and the interval
  
  print *,'Insert the initial time-step.'
  read *, init_ts
  
  print *,'Insert the final time-step.'
  read *, end_ts
  
  print *,'Insert the interval between each time-step.'
  read *, ts_interval
 
  ! Calculate the number of Snapshots
  
  n_snaps = (end_ts - init_ts)/ts_interval 
  
  ! Store the initial time step in the ts counter
  
  ts = init_ts
  
  ! Allocation of memory
  
  allocate(snap(n_of_cells,17));           snap(:,:)=(0.0,0.0)
  allocate(sum_matrix(n_of_cells,15));     sum_matrix(:,:)=(0.0,0.0) 
  allocate(sum_matrix_2(n_of_cells,15));   sum_matrix_2(:,:)=(0.0,0.0)
  
!------------------------------------------------------------------------------! 
!            Dynamic Smagorinsky model Snapshots treatment                     !
!------------------------------------------------------------------------------!
  
  if (answer .eq. 1) then
  
  do i = 1, n_snaps
  
	n_c_tot = 0
  
  	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
  
  	open (unit=1+i,file = trim(filename),form='formatted',status='old')
  	
  	do  
 
  		read(1+i, *,end=10)       snap(n_c_tot+1, 1)  ,   &  ! x-coordinate
  		                          snap(n_c_tot+1, 2)  ,   &  ! y-coordinate
  		                          snap(n_c_tot+1, 3)  ,   &  ! z-coordinate
  		                          snap(n_c_tot+1, 4)  ,   &  ! u-velocity component
  		                          snap(n_c_tot+1, 5)  ,   &  ! v-velocity component
  		                          snap(n_c_tot+1, 6)  ,   &  ! w-velocity component
  		                          snap(n_c_tot+1, 7)  ,   &  ! pressure
  		                          snap(n_c_tot+1, 8)  ,   &  ! du/dx
  		                          snap(n_c_tot+1, 9)  ,   &  ! du/dy
  		                          snap(n_c_tot+1, 10) ,   &  ! du/dz
  		                          snap(n_c_tot+1, 11) ,   &  ! dv/dx
  		                          snap(n_c_tot+1, 12) ,   &  ! dv/dy
  		                          snap(n_c_tot+1, 13) ,   &  ! dv/dz
  		                          snap(n_c_tot+1, 14) ,   &  ! dw/dx
  		                          snap(n_c_tot+1, 15) ,   &  ! dw/dy
  		                          snap(n_c_tot+1, 16) ,   &  ! dw/dz
  		                          snap(n_c_tot+1, 17)        ! ni_sgs  		                          
  		                            		                        		                              		                          
  		sum_matrix(n_c_tot+1, 1)  = sum_matrix(n_c_tot+1, 1)  + snap(n_c_tot+1, 1)  ! x-coordinate
  		sum_matrix(n_c_tot+1, 2)  = sum_matrix(n_c_tot+1, 2)  + snap(n_c_tot+1, 2)  ! y-coordinate
  		sum_matrix(n_c_tot+1, 3)  = sum_matrix(n_c_tot+1, 3)  + snap(n_c_tot+1, 3)  ! z-coordinate
  		sum_matrix(n_c_tot+1, 4)  = sum_matrix(n_c_tot+1, 4)  + snap(n_c_tot+1, 4)  ! U-velocity component
  		sum_matrix(n_c_tot+1, 5)  = sum_matrix(n_c_tot+1, 5)  + snap(n_c_tot+1, 5)  ! V-velocity component
  		sum_matrix(n_c_tot+1, 6)  = sum_matrix(n_c_tot+1, 6)  + snap(n_c_tot+1, 6)  ! W-velocity component
  		sum_matrix(n_c_tot+1, 7)  = sum_matrix(n_c_tot+1, 7)  + snap(n_c_tot+1, 7)  ! avg pressure
  		
  		sij(1,1) = snap(n_c_tot+1, 8) ! S11
  		sij(2,2) = snap(n_c_tot+1,12) ! S22
  		sij(3,3) = snap(n_c_tot+1,16) ! S33
  		
  		sij(1,2) = 0.5*(snap(n_c_tot+1, 9) + snap(n_c_tot+1,11)) ! S12
  		sij(1,3) = 0.5*(snap(n_c_tot+1,10) + snap(n_c_tot+1,14)) ! S13
  		sij(2,3) = 0.5*(snap(n_c_tot+1,13) + snap(n_c_tot+1,15)) ! S23
  		
  		ni_sgs = snap(n_c_tot+1,17)   ! ni_sgs
  		
  		tau(1,1) = ni_sgs*sij(1,1) ! sub-grid stress tensor
  		tau(2,2) = ni_sgs*sij(2,2)
  		tau(3,3) = ni_sgs*sij(3,3)
  		
  		tau(1,2) = ni_sgs*sij(1,2)
  		tau(1,3) = ni_sgs*sij(1,3)
  		tau(2,3) = ni_sgs*sij(2,3)
  		
                sum_matrix(n_c_tot+1, 8)  = sum_matrix(n_c_tot+1, 8)  + tau(1,1)  ! avg sub-grid stress tensor, component 11

                sum_matrix(n_c_tot+1, 9)  = sum_matrix(n_c_tot+1, 9)  + tau(2,2)  ! avg sub-grid stress tensor, component 22

                sum_matrix(n_c_tot+1, 10) = sum_matrix(n_c_tot+1, 10) + tau(3,3)  ! avg sub-grid stress tensor, component 33

                sum_matrix(n_c_tot+1, 11) = sum_matrix(n_c_tot+1, 11) + tau(1,2)  ! avg sub-grid stress tensor, component 12

                sum_matrix(n_c_tot+1, 12) = sum_matrix(n_c_tot+1, 12) + tau(1,3)  ! avg sub-grid stress tensor, component 13

                sum_matrix(n_c_tot+1, 13) = sum_matrix(n_c_tot+1, 13) + tau(2,3)  ! avg sub-grid stress tensor, component 23
                
                sum_matrix(n_c_tot+1, 14) = sum_matrix(n_c_tot+1, 14) + tau(1,1)*sij(1,1) + tau(2,2)*sij(2,2) + tau(3,3)*sij(3,3) &
                                                                   + 2*(tau(1,2)*sij(1,2) + tau(1,3)*sij(1,3) + tau(2,3)*sij(2,3)) 
                                                                   
                                                                   ! sub-grid dissipation (epsilon_sgs)
                                                                   
                sum_matrix(n_c_tot+1, 15) = sum_matrix(n_c_tot+1, 15) + (snap(n_c_tot+1,11) - snap(n_c_tot+1,9)) ! z-vorticity component
 		  		
  		n_c_tot = n_c_tot + 1
  		
  	end do
  	
  	10 close (unit=1+i, status='keep')  ! 'keep' to keep the Snapshots, otherwise 'delete' to delete them
  		
  	ts = ts + ts_interval               ! Update of the ts to open the next Snapshot
  	
  end do

!------------------------------------------------------------------------------! 
!            Tensorial viscosity model Snapshots treatment                     !
!------------------------------------------------------------------------------!                             
  		
  else if (answer .eq. 2) then
  
  	! Inertia file reading
  	
  	allocate(inertia(n_of_cells,10)); inertia(:,:)=(0.0,0.0)
  
  	open (unit=789, file='Inertia.txt', form='formatted', status='old') 
  
  	n_c_tot = 0
  
 	 do 
  
  		read (789, *, end=11) inertia(n_c_tot+1, 1), & ! x-coordinate
  			   	      inertia(n_c_tot+1, 2), & ! y-coordinate
  			    	      inertia(n_c_tot+1, 3), & ! z-coordinate
  			    	      inertia(n_c_tot+1, 4), & ! Ixx
  			    	      inertia(n_c_tot+1, 5), & ! Iyy
  			    	      inertia(n_c_tot+1, 6), & ! Izz
  			    	      inertia(n_c_tot+1, 7), & ! Ixy
  			    	      inertia(n_c_tot+1, 8), & ! Ixz
  			    	      inertia(n_c_tot+1, 9), & ! Iyz
  			    	      inertia(n_c_tot+1,10)    ! cell volume
  			    
  		n_c_tot = n_c_tot + 1
  	
  	end do
  
  	11 close(789, status='keep')
  	
  	! Snapshots reading
  	
  	  do i = 1, n_snaps
    
  	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
  
  	open (unit=1+i,file = trim(filename),form='formatted',status='old')
  	
  	do j = 1, n_c_tot  
 
  		read(1+i, *)              snap(j, 1)  ,   &  ! x-coordinate
  		                          snap(j, 2)  ,   &  ! y-coordinate
  		                          snap(j, 3)  ,   &  ! z-coordinate
  		                          snap(j, 4)  ,   &  ! u-velocity component
  		                          snap(j, 5)  ,   &  ! v-velocity component
  		                          snap(j, 6)  ,   &  ! w-velocity component
  		                          snap(j, 7)  ,   &  ! pressure
  		                          snap(j, 8)  ,   &  ! du/dx
  		                          snap(j, 9)  ,   &  ! du/dy
  		                          snap(j, 10) ,   &  ! du/dz
  		                          snap(j, 11) ,   &  ! dv/dx
  		                          snap(j, 12) ,   &  ! dv/dy
  		                          snap(j, 13) ,   &  ! dv/dz
  		                          snap(j, 14) ,   &  ! dw/dx
  		                          snap(j, 15) ,   &  ! dw/dy
  		                          snap(j, 16)        ! dw/dz
  		                           		                            		                            		                        		                              		                          
  		sum_matrix(j, 1)  = sum_matrix(j, 1)  + snap(j, 1)  ! x-coordinate
  		sum_matrix(j, 2)  = sum_matrix(j, 2)  + snap(j, 2)  ! y-coordinate
  		sum_matrix(j, 3)  = sum_matrix(j, 3)  + snap(j, 3)  ! z-coordinate
  		sum_matrix(j, 4)  = sum_matrix(j, 4)  + snap(j, 4)  ! U-velocity component
  		sum_matrix(j, 5)  = sum_matrix(j, 5)  + snap(j, 5)  ! V-velocity component
  		sum_matrix(j, 6)  = sum_matrix(j, 6)  + snap(j, 6)  ! W-velocity component
  		sum_matrix(j, 7)  = sum_matrix(j, 7)  + snap(j, 7)  ! avg pressure
  		
  		ip(1) = (0.5*(+inertia(j,4) - inertia(j,5) - inertia(j,6)))  ! Ixx'
  		ip(2) = (0.5*(-inertia(j,4) + inertia(j,5) - inertia(j,6)))  ! Iyy'
  		ip(3) = (0.5*(-inertia(j,4) - inertia(j,5) + inertia(j,6)))  ! Izz'
  		ip(4) = inertia(j,7)                                         ! Ixy'
  		ip(5) = inertia(j,8)                                         ! Ixz'
  		ip(6) = inertia(j,9)                                         ! Iyz'
  		
  		cv = inertia(j,10)                                           ! cell volume
  		
  		sgv(1,1) = (0.5/cv)*(ip(1)*snap(j,8)  + ip(4)*snap(j,9)  + ip(5)*snap(j,10)) ! subgrid viscosity tensor, component 11
  		sgv(2,2) = (0.5/cv)*(ip(4)*snap(j,11) + ip(2)*snap(j,12) + ip(6)*snap(j,13)) ! subgrid viscosity tensor, component 22       
  		sgv(3,3) = (0.5/cv)*(ip(5)*snap(j,14) + ip(6)*snap(j,15) + ip(3)*snap(j,16)) ! subgrid viscosity tensor, component 33
  		
  		sgv(1,2) = (0.5/cv)*(ip(1)*snap(j,11) + ip(4)*snap(j,12) + ip(5)*snap(j,13)) ! subgrid viscosity tensor, component 12
  		sgv(2,1) = (0.5/cv)*(ip(4)*snap(j,8)  + ip(2)*snap(j,9)  + ip(6)*snap(j,10)) ! subgrid viscosity tensor, component 21
  		sgv(1,3) = (0.5/cv)*(ip(1)*snap(j,14) + ip(4)*snap(j,15) + ip(5)*snap(j,16)) ! subgrid viscosity tensor, component 13
  		
  		sgv(3,1) = (0.5/cv)*(ip(5)*snap(j,8)  + ip(6)*snap(j,9)  + ip(3)*snap(j,10)) ! subgrid viscosity tensor, component 31
  		sgv(2,3) = (0.5/cv)*(ip(4)*snap(j,14) + ip(2)*snap(j,15) + ip(6)*snap(j,16)) ! subgrid viscosity tensor, component 23
  		sgv(3,2) = (0.5/cv)*(ip(5)*snap(j,11) + ip(6)*snap(j,12) + ip(3)*snap(j,13)) ! subgrid viscosity tensor, component 32
  		
  		tau(1,1) = - 2*sgv(1,1)*snap(j,8)  - 2*sgv(2,1)*snap(j,9)  - 2*sgv(3,1)*snap(j,10)   ! subgrid stress, component 11
  		tau(2,2) = - 2*sgv(1,2)*snap(j,11) - 2*sgv(2,2)*snap(j,12) - 2*sgv(3,2)*snap(j,13)   ! subgrid stress, component 22
  		tau(3,3) = - 2*sgv(1,3)*snap(j,14) - 2*sgv(2,3)*snap(j,15) - 2*sgv(3,3)*snap(j,16)   ! subgrid stress, component 33
  		
  		tau(1,2) = - sgv(1,2)*snap(j,8)  - sgv(1,1)*snap(j,11) &                             ! subgrid stress, component 12
  		           - sgv(2,2)*snap(j,9)  - sgv(2,1)*snap(j,12) &
  		           - sgv(3,2)*snap(j,10) - sgv(3,1)*snap(j,13) 
  		           
  		tau(1,3) = - sgv(1,3)*snap(j,8)  - sgv(1,1)*snap(j,14) &                             ! subgrid stress, component 13
  		           - sgv(2,3)*snap(j,9)  - sgv(2,1)*snap(j,15) &
  		           - sgv(3,3)*snap(j,10) - sgv(3,1)*snap(j,16)
  		           
  		tau(2,3) = - sgv(1,3)*snap(j,11)  - sgv(1,2)*snap(j,14) &                            ! subgrid stress, component 23
  		           - sgv(2,3)*snap(j,12)  - sgv(2,2)*snap(j,15) &
  		           - sgv(3,3)*snap(j,13)  - sgv(3,2)*snap(j,16)
  		           
  		sij(1,1) = snap(j+1, 8) ! S11
  		sij(2,2) = snap(j+1,12) ! S22
  		sij(3,3) = snap(j+1,16) ! S33
  		
  		sij(1,2) = 0.5*(snap(j+1, 9) + snap(j+1,11)) ! S12
  		sij(1,3) = 0.5*(snap(j+1,10) + snap(j+1,14)) ! S13
  		sij(2,3) = 0.5*(snap(j+1,13) + snap(j+1,15)) ! S23
  		
  	        sum_matrix(j+1, 8)  = sum_matrix(j+1, 8)  + tau(1,1)  ! avg sub-grid stress tensor, component 11

                sum_matrix(j+1, 9)  = sum_matrix(j+1, 9)  + tau(2,2)  ! avg sub-grid stress tensor, component 22

                sum_matrix(j+1, 10) = sum_matrix(j+1, 10) + tau(3,3)  ! avg sub-grid stress tensor, component 33

                sum_matrix(j+1, 11) = sum_matrix(j+1, 11) + tau(1,2)  ! avg sub-grid stress tensor, component 12

                sum_matrix(j+1, 12) = sum_matrix(j+1, 12) + tau(1,3)  ! avg sub-grid stress tensor, component 13

                sum_matrix(j+1, 13) = sum_matrix(j+1, 13) + tau(2,3)  ! avg sub-grid stress tensor, component 23
                
                sum_matrix(j+1, 14) = sum_matrix(j+1, 14) + tau(1,1)*sij(1,1) + tau(2,2)*sij(2,2) + tau(3,3)*sij(3,3) &
                                                       + 2*(tau(1,2)*sij(1,2) + tau(1,3)*sij(1,3) + tau(2,3)*sij(2,3))  ! sub-grid dissipation 
                                                                                                                        !(epsilon_sgs)
                                                       
                sum_matrix(n_c_tot+1, 15) = sum_matrix(n_c_tot+1, 15) + (snap(n_c_tot+1,11) - snap(n_c_tot+1,9)) ! z-vorticity component
 		  		  		    		
  		           
   	end do
  	
  	close (unit=1+i, status='keep')	  ! 'keep' to keep the Snapshots, otherwise 'delete' to delete them
  		
  	ts = ts + ts_interval             ! Update of the ts to open the next Snapshot
  	
  end do
  
  deallocate(inertia) 
  
!------------------------------------------------------------------------------! 
!                    Implict LES Snapshots treatment                           !
!------------------------------------------------------------------------------!		           
  		            
  else if (answer .eq. 3) then
  
  do i = 1, n_snaps
  
	n_c_tot = 0
  
  	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
  
  	open (unit=1+i,file = trim(filename),form='formatted',status='old')
  	
  	do  
 
  		read(1+i, *,end=13)       snap(n_c_tot+1, 1)  ,   &  ! x-coordinate
  		                          snap(n_c_tot+1, 2)  ,   &  ! y-coordinate
  		                          snap(n_c_tot+1, 3)  ,   &  ! z-coordinate
  		                          snap(n_c_tot+1, 4)  ,   &  ! u-velocity component
  		                          snap(n_c_tot+1, 5)  ,   &  ! v-velocity component
  		                          snap(n_c_tot+1, 6)  ,   &  ! w-velocity component
  		                          snap(n_c_tot+1, 7)  ,   &  ! pressure
  		                          snap(n_c_tot+1, 8)  ,   &  ! du/dx
  		                          snap(n_c_tot+1, 9)  ,   &  ! du/dy
  		                          snap(n_c_tot+1, 10) ,   &  ! du/dz
  		                          snap(n_c_tot+1, 11) ,   &  ! dv/dx
  		                          snap(n_c_tot+1, 12) ,   &  ! dv/dy
  		                          snap(n_c_tot+1, 13) ,   &  ! dv/dz
  		                          snap(n_c_tot+1, 14) ,   &  ! dw/dx
  		                          snap(n_c_tot+1, 15) ,   &  ! dw/dy
  		                          snap(n_c_tot+1, 16)        ! dw/dz
  		                            		                            		                            		                        		                              		                          
  		sum_matrix(n_c_tot+1, 1)  = sum_matrix(n_c_tot+1, 1)  + snap(n_c_tot+1, 1)  ! x-coordinate
  		sum_matrix(n_c_tot+1, 2)  = sum_matrix(n_c_tot+1, 2)  + snap(n_c_tot+1, 2)  ! y-coordinate
  		sum_matrix(n_c_tot+1, 3)  = sum_matrix(n_c_tot+1, 3)  + snap(n_c_tot+1, 3)  ! z-coordinate
  		sum_matrix(n_c_tot+1, 4)  = sum_matrix(n_c_tot+1, 4)  + snap(n_c_tot+1, 4)  ! U-velocity component
  		sum_matrix(n_c_tot+1, 5)  = sum_matrix(n_c_tot+1, 5)  + snap(n_c_tot+1, 5)  ! V-velocity component
  		sum_matrix(n_c_tot+1, 6)  = sum_matrix(n_c_tot+1, 6)  + snap(n_c_tot+1, 6)  ! W-velocity component
  		sum_matrix(n_c_tot+1, 7)  = sum_matrix(n_c_tot+1, 7)  + snap(n_c_tot+1, 7)  ! avg pressure
  		
  		sum_matrix(n_c_tot+1, 15) = sum_matrix(n_c_tot+1, 15) + (snap(n_c_tot+1,11) - snap(n_c_tot+1,9)) ! z-vorticity component
  		 		  		
  		n_c_tot = n_c_tot + 1
  		
  	end do
  	
  	13 close (unit=1+i, status='keep') ! 'keep' to keep the Snapshots, otherwise 'delete' to delete them
  		
  	ts = ts + ts_interval              ! Update of the ts to open the next Snapshot
  	
  end do
		                                            		 		
  end if
  
!------------------------------------------------------------------------------!

  ! Calculate the reduced number of cells due to the spanwise average
  
  n_c_tot_red = n_c_tot/n_layers ! Always possible being n_c_tot a multiple of n_layers (extruded mesh).
  
  ! Calculation of the spanwise average
 
  do j=0, n_c_tot, n_layers
  	
  	do i=1, n_layers
  	
  		sum_matrix_2(a,1)  = sum_matrix_2(a,1)  + sum_matrix(j+i,1)  ! x-coordinate
  		sum_matrix_2(a,2)  = sum_matrix_2(a,2)  + sum_matrix(j+i,2)  ! y-coordinate
  		sum_matrix_2(a,3)  = 0.5*n_layers                            ! If we skip to add the time_average matrix, 
  		                                                             !the z-coordinate is 0.5 as default 
  		                                                             !(useful as a standard for post_processing)
  		                                                             
  		sum_matrix_2(a,4)  = sum_matrix_2(a,4)  + sum_matrix(j+i,4)  ! U velocity component
  		sum_matrix_2(a,5)  = sum_matrix_2(a,5)  + sum_matrix(j+i,5)  ! V-velocity component
  		sum_matrix_2(a,6)  = sum_matrix_2(a,6)  + sum_matrix(j+i,6)  ! W-velocity component
  		sum_matrix_2(a,7)  = sum_matrix_2(a,7)  + sum_matrix(j+i,7)  ! average pressure P
  		
  		sum_matrix_2(a,15) = sum_matrix_2(a,15) + sum_matrix(j+i,15) ! z-vorticity component
  		
  	if (answer .eq. 1 .or. answer .eq. 2) then	
  		
  		sum_matrix_2(a,8)  = sum_matrix_2(a,8)  + sum_matrix(j+i,8)  ! sub-grid stress tensor, component 11 
  		sum_matrix_2(a,9)  = sum_matrix_2(a,9)  + sum_matrix(j+i,9)  ! sub-grid stress tensor, component 22
  		sum_matrix_2(a,10) = sum_matrix_2(a,10) + sum_matrix(j+i,10) ! sub-grid stress tensor, component 33
  		sum_matrix_2(a,11) = sum_matrix_2(a,11) + sum_matrix(j+i,11) ! sub-grid stress tensor, component 12
  		sum_matrix_2(a,12) = sum_matrix_2(a,12) + sum_matrix(j+i,12) ! sub-grid stress tensor, component 13
  		sum_matrix_2(a,13) = sum_matrix_2(a,13) + sum_matrix(j+i,13) ! sub-grid stress tensor, component 23
  		sum_matrix_2(a,14) = sum_matrix_2(a,14) + sum_matrix(j+i,14) ! sub-grid dissipation (epsilon_sgs)
  		  		  		  		  		  		   		  		
  	end if 		
  		 			
  	end do
  	
  	a = a + 1      ! This parameter allows to store the results in the summation matrix in a reduced number of rows.
                       ! Final total 'a' will be equal to n_c_tot_red
  end do
  
!------------------------------------------------------------------------------!
  
  ! Dividing the sum matrix by the number time steps and the number of layers (averaging)
  
  sum_matrix_2 = sum_matrix_2/n_snaps
    
  sum_matrix_2 = sum_matrix_2/n_layers
  	    
!------------------------------------------------------------------------------!
  
  ! Calculation of the Reynolds stresses, <u'u'>, <v'v'>, <w',w'>, <u'v'> 
    
  allocate(u_prime(n_c_tot,n_snaps));  u_prime(:,:)=(0.0,0.0)
  allocate(v_prime(n_c_tot,n_snaps));  v_prime(:,:)=(0.0,0.0)
  allocate(w_prime(n_c_tot,n_snaps));  w_prime(:,:)=(0.0,0.0)
  
  allocate(uu(n_c_tot,n_snaps));       uu(:,:)=(0.0,0.0)
  allocate(vv(n_c_tot,n_snaps));       vv(:,:)=(0.0,0.0)
  allocate(ww(n_c_tot,n_snaps));       ww(:,:)=(0.0,0.0)
  allocate(rs(n_c_tot,n_snaps));       rs(:,:)=(0.0,0.0)
  
  allocate(sum_uu(n_c_tot));           sum_uu(:)=(0.0)
  allocate(sum_vv(n_c_tot));           sum_vv(:)=(0.0)
  allocate(sum_ww(n_c_tot));           sum_ww(:)=(0.0)
  allocate(sum_rs(n_c_tot));           sum_rs(:)=(0.0)
  
  ! Re-set the ts to the initial one and the parameter a to 1
  
  ts = init_ts
  
  a = 1
  
  ! DSM
  
  if (answer .eq. 1) then
    
  do i = 1, n_snaps
  
  	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
  
  	open (unit=1+i,file = trim(filename),form='formatted',status='old')
  	
  	do j = 1, n_c_tot 
  		
  		read(1+i, *)              snap(j, 1)  ,   &  ! x-coordinate
  		                          snap(j, 2)  ,   &  ! y-coordinate
  		                          snap(j, 3)  ,   &  ! z-coordinate
  		                          snap(j, 4)  ,   &  ! u-velocity component
  		                          snap(j, 5)  ,   &  ! v-velocity component
  		                          snap(j, 6)  ,   &  ! w-velocity component
  		                          snap(j, 7)  ,   &  ! p, pressure
  		                          snap(j, 8)  ,   &  ! du/dx
  		                          snap(j, 9)  ,   &  ! du/dy
  		                          snap(j, 10) ,   &  ! du/dz
  		                          snap(j, 11) ,   &  ! dv/dx
  		                          snap(j, 12) ,   &  ! dv/dy
  		                          snap(j, 13) ,   &  ! dv/dz
  		                          snap(j, 14) ,   &  ! dw/dx
  		                          snap(j, 15) ,   &  ! dw/dy
  		                          snap(j, 16) ,   &  ! dw/dz
  		                          snap(j, 17)        ! ni_sgs
  		                           		               		             
  		u_prime(j,i) = snap(j,4) - sum_matrix_2(a,4)   ! u' fluctuations for each row
  		v_prime(j,i) = snap(j,5) - sum_matrix_2(a,5)   ! v' fluctuations for each row
  		w_prime(j,i) = snap(j,6) - sum_matrix_2(a,6)   ! w' fluctuations for each row
  		
  		uu(j,i) = u_prime(j,i)*u_prime(j,i)          ! Calculation of the instanteneous Reynolds stress u'u'
  		vv(j,i) = v_prime(j,i)*v_prime(j,i)          ! Calculation of the instanteneous Reynolds stress v'v'
  		ww(j,i) = w_prime(j,i)*w_prime(j,i)          ! Calculation of the instanteneous Reynolds stress w'w'
  		rs(j,i) = u_prime(j,i)*v_prime(j,i)          ! Calculation of the instanteneous Reynolds stress u'v'
  		
  		sum_uu(j) = sum_uu(j) + uu(j,i) 
  		sum_vv(j) = sum_vv(j) + vv(j,i) 
  		sum_ww(j) = sum_ww(j) + ww(j,i)              	                          
  		sum_rs(j) = sum_rs(j) + rs(j,i)
  		
  		if (j .eq. a*n_layers) then
  		
  		a = a + 1
  		
  		end if	
  		                                    
        end do
        
        close (unit=1+i, status='keep')
        
        ts = ts + ts_interval
  
  end do
  
  ! TVM or ILES
  
  else if(answer .eq. 2 .or. answer .eq. 3) then
  
  do i = 1, n_snaps
  
  	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
  
  	open (unit=1+i,file = trim(filename),form='formatted',status='old')
  	
  	do j = 1, n_c_tot 
  		
  		read(1+i, *)              snap(j, 1)  ,   &  ! x-coordinate
  		                          snap(j, 2)  ,   &  ! y-coordinate
  		                          snap(j, 3)  ,   &  ! z-coordinate
  		                          snap(j, 4)  ,   &  ! u-velocity component
  		                          snap(j, 5)  ,   &  ! v-velocity component
  		                          snap(j, 6)  ,   &  ! w-velocity component
  		                          snap(j, 7)  ,   &  ! p, pressure
  		                          snap(j, 8)  ,   &  ! du/dx
  		                          snap(j, 9)  ,   &  ! du/dy
  		                          snap(j, 10) ,   &  ! du/dz
  		                          snap(j, 11) ,   &  ! dv/dx
  		                          snap(j, 12) ,   &  ! dv/dy
  		                          snap(j, 13) ,   &  ! dv/dz
  		                          snap(j, 14) ,   &  ! dw/dx
  		                          snap(j, 15) ,   &  ! dw/dy
  		                          snap(j, 16)        ! dw/dz  		                         
  		               		             
  		u_prime(j,i) = snap(j,4) - sum_matrix_2(a,4)   ! u' fluctuations for each row
  		v_prime(j,i) = snap(j,5) - sum_matrix_2(a,5)   ! v' fluctuations for each row
  		w_prime(j,i) = snap(j,6) - sum_matrix_2(a,6)   ! w' fluctuations for each row
  		
  		uu(j,i) = u_prime(j,i)*u_prime(j,i)          ! Calculation of the instanteneous Reynolds stress u'u'
  		vv(j,i) = v_prime(j,i)*v_prime(j,i)          ! Calculation of the instanteneous Reynolds stress v'v'
  		ww(j,i) = w_prime(j,i)*w_prime(j,i)          ! Calculation of the instanteneous Reynolds stress w'w'
  		rs(j,i) = u_prime(j,i)*v_prime(j,i)          ! Calculation of the instanteneous Reynolds stress u'v'
  		
  		sum_uu(a) = sum_uu(a) + uu(j,i) 
  		sum_vv(a) = sum_vv(a) + vv(j,i) 
  		sum_ww(a) = sum_ww(a) + ww(j,i)              	                          
  		sum_rs(a) = sum_rs(a) + rs(j,i)
  		
  		if (j .eq. a*n_layers) then
  		
  		a = a + 1
  		
  		end if                           
        
        end do
        
        close (unit=1+i, status='keep')
        
        ts = ts + ts_interval
  
  end do
  
  end if
  
  ! Time-average and spanwise-average of the Reynolds stress

  sum_uu = sum_uu/n_snaps
  sum_vv = sum_vv/n_snaps
  sum_ww = sum_ww/n_snaps  
  sum_rs = sum_rs/n_snaps
  
  sum_uu = sum_uu/n_layers
  sum_vv = sum_vv/n_layers
  sum_ww = sum_ww/n_layers  
  sum_rs = sum_rs/n_layers
   
  ! Calculating the turbulent kinetic energy (TKE)
  
  allocate(TKE_vector(n_c_tot_red)); TKE_vector(:)=0.0
  
  j = 0
  
  do j = 1, n_c_tot_red
  
  	TKE_vector(j) = TKE_vector(j) + sum_uu(j) + sum_vv(j) + sum_ww(j)
  	
  end do
                        	
!------------------------------------------------------------------------------! 
  
  ! Writing of the output file 'Time+Span_Average.txt'                                        
  
  ! DSM or TVM (subgrid stress and subgrid dissipation can be displayed)
  
  if (answer .eq. 1 .or. answer .eq. 2) then
    
  open (unit=2, file = 'Time+Span_Average.txt', form='formatted', status='unknown')
  
  	do j = 1, n_c_tot_red + 1
  	
  		if (j .eq. 1) then
  		
  		write (2, *) 'x', ',', 'y', ',', 'z', ',', &
  			     'U', ',', 'V', ',', 'W', ',', &
  			     'P', ',', &
  			     '<u"u">', ',', '<v"v">', ',', &
  			     '<w"w">', ',', '<u"v">', ',', &
  			     '<k>', ',', &
  			     '<tau_11>', ',', '<tau_22>', ',', '<tau_33>', ',', &
  			     '<tau_12>', ',', '<tau_13>', ',', '<tau_23>', ',', &
  			     '<eps_sgs>', ',', '<omega_z>'
  		
  		else
    
  		write (2, *)    sum_matrix_2(j-1, 1), ',', & ! x-coordinate     
  	             		sum_matrix_2(j-1, 2), ',', & ! y-coordinate
  	             		sum_matrix_2(j-1, 3), ',', & ! z-coordinate
  	             		sum_matrix_2(j-1, 4), ',', & ! U-velocity component
  	             		sum_matrix_2(j-1, 5), ',', & ! V-velocity component
  	             		sum_matrix_2(j-1, 6), ',', & ! W-velocity component
  	             		sum_matrix_2(j-1, 7), ',', & ! P, mean pressure
  	             		sum_uu(j-1), ',', &          ! <u"u">
  	             		sum_vv(j-1), ',', &          ! <v"v">
  	             		sum_ww(j-1), ',', &          ! <w"w">	             
  	             		sum_rs(j-1), ',', &          ! <u"v">
  	             		TKE_vector(j-1), ',', &      ! <k>
                   	        sum_matrix_2(j-1, 8), ',', & ! mean subgrid stress tensor, component 11
                                sum_matrix_2(j-1, 9), ',', & ! mean subgrid stress tensor, component 22
                                sum_matrix_2(j-1,10), ',', & ! mean subgrid stress tensor, component 33
                                sum_matrix_2(j-1,11), ',', & ! mean subgrid stress tensor, component 12
                                sum_matrix_2(j-1,12), ',', & ! mean subgrid stress tensor, component 13
                                sum_matrix_2(j-1,13), ',', & ! mean subgrid stress tensor, component 23
                                sum_matrix_2(j-1,14), ',', & ! mean subgrid dissipation (epsilon_sgs)
                                sum_matrix_2(j-1,15)         ! mean z-vorticity component
                                
                end if                      
                                                                             
  	end do
  	
  	close(2)
  	
  ! ILES
  	             		
  else if (answer .eq. 3) then
  
    open (unit=2, file = 'Time+Span_Average.txt', form='formatted', status='unknown')
  
  	do j = 1, n_c_tot_red + 1
  	
  	if (j .eq. 1) then
  		
  		write (2, *) 'x', ',', 'y', ',', 'z', ',', &
  			     'U', ',', 'V', ',', 'W', ',', &
  			     'P', ',', &
  			     '<u"u">', ',', '<v"v">', ',', &
  			     '<w"w">', ',', '<u"v">', ',', &
  			     '<k>', ',', &
  			     '<omega_z>'
  		
  		else
    
  		write (2, *)  sum_matrix_2(j-1, 1), ',', & ! x-coordinate     
  	             		sum_matrix_2(j-1, 2), ',', & ! y-coordinate
  	             		sum_matrix_2(j-1, 3), ',', & ! z-coordinate
  	             		sum_matrix_2(j-1, 4), ',', & ! U-velocity component
  	             		sum_matrix_2(j-1, 5), ',', & ! V-velocity component
  	             		sum_matrix_2(j-1, 6), ',', & ! W-velocity component
  	             		sum_matrix_2(j-1, 7), ',', & ! P, mean pressure
  	             		sum_uu(j-1), ',', &          ! <u"u">
  	             		sum_vv(j-1), ',', &          ! <v"v">
  	             		sum_ww(j-1), ',', &          ! <w"w">	             
  	             		sum_rs(j-1), ',', &          ! <u"v">
  	             		TKE_vector(j-1), ',', &      ! <k>
  	             		sum_matrix_2(j-1,15)         ! mean z-vorticity component
  	             		
  	end if

  	end do
  	
  	close(2)
   	        
  end if
              		 	             		  	                                	           
  deallocate(snap)
  deallocate(sum_matrix)
  deallocate(sum_matrix_2)
  deallocate(TKE_vector)
  deallocate(u_prime)
  deallocate(v_prime)
  deallocate(w_prime)
  deallocate(uu)
  deallocate(vv)
  deallocate(ww)
  deallocate(rs)
  deallocate(sum_uu)
  deallocate(sum_vv)
  deallocate(sum_ww)
  deallocate(sum_rs)
  
!------------------------------------------------------------------------------!
  
  else if (selection .eq. 2) then ! Only spanwise average of a specific Snapshot
  
  print *,'Insert the Snapshot corresponding ts'
  read *, ts
  
  ! Allocation of memory
  
  allocate(snap(n_of_cells,17));         snap(:,:)=(0.0,0.0)
  allocate(sum_matrix(n_of_cells,8));    sum_matrix(:,:)=(0.0,0.0)
  
  ! Creation of the Snapshot filename
  
    	select case(ts)
  	
  	case (0:9)
        write (filename, "(A16,I1,A4)") "SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A15,I2,A4)") "SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A14,I3,A4)") "SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A13,I4,A4)") "SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A12,I5,A4)") "SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A11,I6,A4)") "SnapShot-ts", ts, '.txt'
    end select
 
  ! Snapshot opening
  
   if (answer .eq. 1) then ! Dynamic Smagorinsky 
   
  	open (unit=1,file = trim(filename),form='formatted',status='old')
  	
  	do  
 
  		read(1, *,end=18)         snap(n_c_tot+1, 1)  ,   &  ! x-coordinate
  		                          snap(n_c_tot+1, 2)  ,   &  ! y-coordinate
  		                          snap(n_c_tot+1, 3)  ,   &  ! z-coordinate
  		                          snap(n_c_tot+1, 4)  ,   &  ! u-velocity component
  		                          snap(n_c_tot+1, 5)  ,   &  ! v-velocity component
  		                          snap(n_c_tot+1, 6)  ,   &  ! w-velocity component
  		                          snap(n_c_tot+1, 7)  ,   &  ! pressure
  		                          snap(n_c_tot+1, 8)  ,   &  ! du/dx
  		                          snap(n_c_tot+1, 9)  ,   &  ! du/dy
  		                          snap(n_c_tot+1, 10) ,   &  ! du/dz
  		                          snap(n_c_tot+1, 11) ,   &  ! dv/dx
  		                          snap(n_c_tot+1, 12) ,   &  ! dv/dy
  		                          snap(n_c_tot+1, 13) ,   &  ! dv/dz
  		                          snap(n_c_tot+1, 14) ,   &  ! dw/dx
  		                          snap(n_c_tot+1, 15) ,   &  ! dw/dy
  		                          snap(n_c_tot+1, 16) ,   &  ! dw/dz
  		                          snap(n_c_tot+1, 17)        ! ni_sgs 
  		                          
	n_c_tot = n_c_tot + 1
  		
  	end do
  	
  	18 close (unit=1, status='keep')  ! 'keep' to keep the Snapshot, otherwise 'delete' to delete it
  	
  else if (answer .eq. 2 .or. answer .eq. 3) then ! Tensorial viscosity model or ILES
  
    	open (unit=1,file = trim(filename),form='formatted',status='old')
  	
  	do  
 
  		read(1, *,end=12)         snap(n_c_tot+1, 1)  ,   &  ! x-coordinate
  		                          snap(n_c_tot+1, 2)  ,   &  ! y-coordinate
  		                          snap(n_c_tot+1, 3)  ,   &  ! z-coordinate
  		                          snap(n_c_tot+1, 4)  ,   &  ! u-velocity component
  		                          snap(n_c_tot+1, 5)  ,   &  ! v-velocity component
  		                          snap(n_c_tot+1, 6)  ,   &  ! w-velocity component
  		                          snap(n_c_tot+1, 7)  ,   &  ! pressure
  		                          snap(n_c_tot+1, 8)  ,   &  ! du/dx
  		                          snap(n_c_tot+1, 9)  ,   &  ! du/dy
  		                          snap(n_c_tot+1, 10) ,   &  ! du/dz
  		                          snap(n_c_tot+1, 11) ,   &  ! dv/dx
  		                          snap(n_c_tot+1, 12) ,   &  ! dv/dy
  		                          snap(n_c_tot+1, 13) ,   &  ! dv/dz
  		                          snap(n_c_tot+1, 14) ,   &  ! dw/dx
  		                          snap(n_c_tot+1, 15) ,   &  ! dw/dy
  		                          snap(n_c_tot+1, 16)        ! dw/dz
  		                            		                          
	n_c_tot = n_c_tot + 1
  		
  	end do
  	
  	12 close (unit=1, status='keep')  ! 'keep' to keep the Snapshot, otherwise 'delete' to delete it
  	
  end if
  
  ! Calculate the reduced number of cells due to the spanwise average
  
  n_c_tot_red = n_c_tot/n_layers ! Always possible being n_c_tot a multiple of n_layers (extruded mesh).
  
  ! Calculation of the spanwise average
 
  do j=0, n_c_tot, n_layers
  	
  	do i=1, n_layers
  	
  		sum_matrix(a,1)  = sum_matrix(a,1)  + snap(j+i,1)         ! x-coordinate
  		sum_matrix(a,2)  = sum_matrix(a,2)  + snap(j+i,2)         ! y-coordinate
  		sum_matrix(a,3)  = 0.5*n_layers                           ! z-coordinate (set as reference 0.5)
  		sum_matrix(a,4)  = sum_matrix(a,4)  + snap(j+i,4)         ! U velocity component
  		sum_matrix(a,5)  = sum_matrix(a,5)  + snap(j+i,5)         ! V-velocity component
  		sum_matrix(a,6)  = sum_matrix(a,6)  + snap(j+i,6)         ! W-velocity component
  		sum_matrix(a,7)  = sum_matrix(a,7)  + snap(j+i,7)         ! P, average pressure
   		sum_matrix(a,8)  = sum_matrix(a,8)  + &
  				 + (snap(j+1,11) - snap(j+1,9))	          ! z-vorticity component	
  		 			 
  	end do
  	
  	a = a + 1 ! This parameter allows to store the results in the summation matrix in a reduced number of rows.
  
  end do
  
  ! Dividing the sum matrix by the number of layers (averaging)
  
  sum_matrix = sum_matrix/n_layers
  
  ! Calculating the "spatial Reynolds stresses"
    
  allocate(u_prime(n_c_tot,1));  u_prime(:,:)=(0.0,0.0)
  allocate(v_prime(n_c_tot,1));  v_prime(:,:)=(0.0,0.0)
  allocate(w_prime(n_c_tot,1));  w_prime(:,:)=(0.0,0.0)
  
  allocate(uu(n_c_tot,1));       uu(:,:)=(0.0,0.0)
  allocate(vv(n_c_tot,1));       vv(:,:)=(0.0,0.0)
  allocate(ww(n_c_tot,1));       ww(:,:)=(0.0,0.0)
  allocate(rs(n_c_tot,1));       rs(:,:)=(0.0,0.0)
  
  allocate(sum_uu(n_c_tot));     sum_uu(:)=(0.0)
  allocate(sum_vv(n_c_tot));     sum_vv(:)=(0.0)
  allocate(sum_ww(n_c_tot));     sum_ww(:)=(0.0)
  allocate(sum_rs(n_c_tot));     sum_rs(:)=(0.0)
  
  ! Calculations
  
  a = 1
  
   	do j = 1, n_c_tot 
  		  		                           		               		             
  		u_prime(j,1) = snap(j,4) - sum_matrix(a,4)   ! u' fluctuations for each row
  		v_prime(j,1) = snap(j,5) - sum_matrix(a,5)   ! v' fluctuations for each row
  		w_prime(j,1) = snap(j,6) - sum_matrix(a,6)   ! w' fluctuations for each row
  		
  		uu(j,1) = u_prime(j,1)*u_prime(j,1)          ! Calculation of the instanteneous "Reynolds stress" u'u'
  		vv(j,1) = v_prime(j,1)*v_prime(j,1)          ! Calculation of the instanteneous "Reynolds stress" v'v'
  		ww(j,1) = w_prime(j,1)*w_prime(j,1)          ! Calculation of the instanteneous "Reynolds stress" w'w'
  		rs(j,1) = u_prime(j,1)*v_prime(j,1)          ! Calculation of the instanteneous "Reynolds stress" u'v'
  		
  		sum_uu(j) = sum_uu(j) + uu(j,1) 
  		sum_vv(j) = sum_vv(j) + vv(j,1) 
  		sum_ww(j) = sum_ww(j) + ww(j,1)              	                          
  		sum_rs(j) = sum_rs(j) + rs(j,1)
  		
  		if (j .eq. a*n_layers) then
  		
  		a = a + 1
  		
  		end if	
  		                                    
        end do
         
  ! Spanwise-average of the "Reynolds stresses"
  
  sum_uu = sum_uu/n_layers
  sum_vv = sum_vv/n_layers
  sum_ww = sum_ww/n_layers  
  sum_rs = sum_rs/n_layers
   
  ! Calculating the turbulent kinetic energy (TKE)
  
  allocate(TKE_vector(n_c_tot_red)); TKE_vector(:)=0.0
  
  j = 0
  
  do j = 1, n_c_tot_red
  
  	TKE_vector(j) = TKE_vector(j) + sum_uu(j) + sum_vv(j) + sum_ww(j)
  	
  end do
    
  ! Writing of the output file 'Span_Avg_Snapshot-ts*.txt'
  
  j = 0
  
    ! Creation of the Snapshot z-average filename
  
    	select case(ts)
  	
  	case (0:9)
        write (filename, "(A25,I1,A4)") "Span_Avg_SnapShot-ts00000", ts, '.txt'
      case (10:99)
        write (filename, "(A24,I2,A4)") "Span_Avg_SnapShot-ts0000", ts, '.txt'
      case (100:999)
        write (filename, "(A23,I3,A4)") "Span_Avg_SnapShot-ts000", ts, '.txt'
      case (1000:9999)
        write (filename, "(A22,I4,A4)") "Span_Avg_SnapShot-ts00", ts, '.txt'
      case (10000:99999)
        write (filename, "(A21,I5,A4)") "Span_Avg_SnapShot-ts0", ts, '.txt'
      case default
        write (filename, "(A22,I6,A4)") "Span_Avg_SnapShot-ts", ts, '.txt'
    end select
    
  open (unit=2, file = trim(filename), form='formatted', status='unknown')
  
  do j = 1, n_c_tot_red + 1
    
  		if (j .eq. 1) then
  			     	
  			     	write (2, *) 'x', ',', 'y', ',', 'z', ',', &
  			     		     'U', ',', 'V', ',', 'W', ',', &
  			     	             'P', ',', '<omega_z>', ',', &
  			     	             '<u"u">|z', ',', '<v"v">|z', ',', &
  			                     '<w"w">|z', ',', '<u"v">|z', ',', &
  			                     '<k>|z'
  			     	            			     	  	
  		else	
          		
          		write (2, *) sum_matrix(j-1, 1), ',', &
          		     	     sum_matrix(j-1, 2), ',', &
          		     	     sum_matrix(j-1, 3), ',', &
          		             sum_matrix(j-1, 4), ',', &
          		             sum_matrix(j-1, 5), ',', &
          		             sum_matrix(j-1, 6), ',', &
          		             sum_matrix(j-1, 7), ',', &
          		             sum_matrix(j-1, 8), ',', &
          		             sum_uu(j-1), ',', &
          		             sum_vv(j-1), ',', &
          		             sum_ww(j-1), ',', &
          		             sum_rs(j-1), ',', &
          		             TKE_vector(j-1)
          		          		         		          	         		             
                end if
   	             
  end do 
  
  deallocate(sum_matrix)
  deallocate(snap)
  deallocate(TKE_vector)
  deallocate(u_prime)
  deallocate(v_prime)
  deallocate(w_prime)
  deallocate(uu)
  deallocate(vv)
  deallocate(ww)
  deallocate(rs)
  deallocate(sum_uu)
  deallocate(sum_vv)
  deallocate(sum_ww)
  deallocate(sum_rs)
    
  close(2)
      		                            		                            		                        		                              		                       
  end if ! Closure of the main if
     
!------------------------------------------------------------------------------!
 
  print *,'Done!'
 
  end program  
  
!------------------------------------------------------------------------------!
  
  subroutine logo()
  
  print *,'#===================================' // &
          '===================================='
  print *,'#'
  print *,'#                     ___________________  ___      '
  print *,'#                     |  ___//\__    __// / _ \\    '
  print *,'#                     |  |       |  ||   / /_\ \\   '
  print *,'#                     |  |__     |  ||  /  ___  \\  '
  print *,'#                     |_____\\   |__|| /__/   \__\\ '
  print *,'#                                                   '
  print *,'#                      Calculate_Total_Average.f95   '
  print *,'#                                                   '
  print *,'#         This program calculates the total-average of the Snapshots,'
  print *,'#               so the time and the spanwise (z) averages.'
  print *,'#'
  print *,'#      If selected, this program can handle also the spanwise average'
  print *,'#                       of a specific Snapshot.'
  print *,'#'
  print *,'#           The output file is ordered in the following manner:'
  print *,'#               x,y,z,U,V,W,P, <u"u">, <v"v">, <w"w">, <u"v">, '
  print *,'#        <tau_11>, <tau_22>, <tau_33>, <tau_12>, <tau_13>, <tau_23>,'
  print *,'#                         <eps_sgs>, <omega_z>'
  print *,'#'
  print *,'#        Subgrid stress and dissipation are not written for ILES.'
  print *,'#'
  print *,'#         For spanwise-average only, a sort of Reynolds stresses' 
  print *,'#                          are calculated.'
  print *,'#'
  print *,'#'                                                  
  print *,'#===================================' // &
          '===================================='
  
  return
  end
  
  
  
  
