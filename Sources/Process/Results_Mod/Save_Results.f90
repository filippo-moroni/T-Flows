!==============================================================================!
  subroutine Save_Results(Results,  &
                          Flow, Turb, Vof, Swarm, ts, plot_inside, domain)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
!   and .txt format (Snapshots)                                                !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)       :: Results
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  integer                   :: ts           ! time step
  logical                   :: plot_inside  ! plot results inside?
  integer,         optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Var_Type),  pointer     :: phi
  integer(SP)                  :: data_size
  integer                      :: data_offset, cell_offset, i_fac
  integer                      :: s, n, n_conns, n_polyg, sc, f8, f9, ua, run
  integer                      :: s1, s2, c1, c2, c_f, c_l
  real                         :: dist1, dist2
  character(SL)                :: name_out_8, name_out_9, name_mean, a_name
  character(SL)                :: str1, str2
  integer, pointer, contiguous :: int_save(:), type_save(:), offs_save(:)
  real,    pointer, contiguous :: save_01(:), save_02(:), save_03(:)
  real,    pointer, contiguous :: save_04(:), save_05(:), save_06(:)
  real,    pointer, contiguous :: var_ins(:)
  real,    pointer, contiguous :: v2_calc(:), kin_vis_t(:), phi_save(:)
  
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

!------------------------------------------------------------------------------!
!                        Snapshots area by Yves (1/2)                          !
!------------------------------------------------------------------------------!
  
  type SnapShot	

    real, allocatable :: u_snap(:), v_snap(:), w_snap(:), p_snap(:)	! x,y and z velocity components and pressure	  
    real, allocatable :: x_snap(:), y_snap(:), z_snap(:)		! x,y and z coordinates of the cell
    real, allocatable :: ni_turb(:)                                     ! generic turbulent viscosity	
    
    real, allocatable :: dudx_snap(:), dudy_snap(:), dudz_snap(:)       ! gradient components for the u velocity component
    real, allocatable :: dvdx_snap(:), dvdy_snap(:), dvdz_snap(:)	! gradient components for the v velocity component
    real, allocatable :: dwdx_snap(:), dwdy_snap(:), dwdz_snap(:)       ! gradient components for the w velocity component
    
    integer, allocatable :: pos_ind(:)     ! Position index of the cell in the ordered single vector
    
    real,    allocatable :: x_2(:)         ! "x" coordinate of the cell. It is necessary, otherwise "z_snap" would be ordered differently after sorting
    real,    allocatable :: y_2(:)	   ! "y" coordinate of the cell. ""  ""
    real,    allocatable :: z_2(:)	   ! "z" coordinate of the cell. ""  ""
    
  end type
  
  type(SnapShot) :: Snap		! Creation of the variable "Snap", of type "SnapShot"
  
  ! Local variables definition
  
  integer :: i = 0		
  integer :: j = 0			
  integer :: n_of_cells = 10000000	! This integer must be larger than the number of cells of the domain; it is used to allocate the memory
  
  character(len=1024) :: filename	! Filename of the file to open or to close
  
  integer :: n_c_tot = 0  		! Counter for the total number of cells

!------------------------------------------------------------------------------!

  call Profiler % Start('Save_Vtu_Results')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Take aliases
  Grid => Flow % pnt_grid

  if(.not. plot_inside .and. .not. Results % boundary) return

  call Work % Connect_Int_Cell(int_save, type_save, offs_save)
  call Work % Connect_Real_Cell(save_01, save_02, save_03,  &
                                save_04, save_05, save_06)
  call Work % Connect_Real_Cell(var_ins, v2_calc, kin_vis_t, phi_save)

  !------------------------------------------------!
  !   Mark the beginnings and end of cell ranges   !
  !------------------------------------------------!

  ! For cells inside
  if(plot_inside) then
    c_f = 1
    c_l = Grid % n_cells
    if(.not. PLOT_BUFFERS) c_l = Grid % n_cells - Grid % Comm % n_buff_cells

  ! For boundary cells
  else
    c_f = -Grid % n_bnd_cells
    c_l = -1
    if(.not. PLOT_BUFFERS) then
      do c_f = -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_f) .eq. this_proc) exit
      end do
      do c_l = -1, -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_l) .eq. this_proc) exit
      end do
    end if
  end if

  !-------------------------------------------------------------------------!
  !   Count connections and polygons in this Grid, you will need it later   !
  !-------------------------------------------------------------------------!
  n_conns = 0
  n_polyg = 0
  if(plot_inside) then
    ! Connections
    do c1 = c_f, c_l
      n_conns = n_conns + abs(Grid % cells_n_nodes(c1))
    end do
    ! Polygons
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                   ! add one for number of polyfs
        do i_fac = 1, Grid % cells_n_faces(c1)  ! add all faces and their nodes
          s = Grid % cells_f(i_fac, c1)
          n = Grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end do
  else
    do c2 = c_f, c_l
      n_conns = n_conns + Grid % cells_n_nodes(c2)
    end do
  end if

  call Comm_Mod_Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File % Set_Name(name_out_8,             &
                         time_step=ts,           &
                         extension='.pvtu',      &
                         domain=domain)
    call File % Set_Name(name_out_9,             &
                         processor=this_proc,    &
                         time_step=ts,           &
                         extension='.vtu',       &
                         domain=domain)
  else
    call File % Set_Name(name_out_8,             &
                         time_step=ts,           &
                         appendix ='-bnd',       &
                         extension='.pvtu',      &
                         domain=domain)
    call File % Set_Name(name_out_9,             &
                         processor=this_proc,    &
                         time_step=ts,           &
                         appendix ='-bnd',       &
                         extension='.vtu',       &
                         domain=domain)
  end if

  if(n_proc > 1 .and. this_proc .eq. 1) then
    call File % Open_For_Writing_Binary(name_out_8, f8)
  end if
  call File % Open_For_Writing_Binary(name_out_9, f9)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_0 // '<?xml version="1.0"?>'              // LF
    write(f8) IN_0 // '<VTKFile type="PUnstructuredGrid">' // LF
    write(f8) IN_1 // '<PUnstructuredGrid GhostLevel="1">' // LF
  end if

  write(f9) IN_0 // '<?xml version="1.0"?>'                           // LF
  write(f9) IN_0 // '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                    'byte_order="LittleEndian">'                      // LF
  write(f9) IN_1 // '<UnstructuredGrid>'                              // LF

  write(str1,'(i0.0)') Grid % n_nodes
  if(plot_inside) then
    write(str2,'(i0.0)') (c_l-c_f+1)
  else
    if((c_l-c_f+1) .eq. 0) then
      write(str2,'(i1)')   (c_l-c_f+1)  ! 0.0 doesn't work for zero :-/
    else
      write(str2,'(i0.0)') (c_l-c_f+1)
    end if
  end if
  write(f9) IN_2 // '<Piece NumberOfPoints="' // trim(str1)    //  &
                    '" NumberOfCells ="' // trim(str2) // '">' // LF

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_3 // '<PPoints>' // LF
    write(f8) IN_4 // '<PDataArray type='//floatp  //  &
                      ' NumberOfComponents="3"/>'  // LF
    write(f8) IN_3 // '</PPoints>' // LF
  end if

  write(str1, '(i1)') data_offset
  write(f9) IN_3 // '<Points>'                        // LF
  write(f9) IN_4 // '<DataArray type='//floatp        //  &
                    ' NumberOfComponents="3"'         //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
  write(f9) IN_4 // '</DataArray>' // LF
  write(f9) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3

  !-----------!
  !   Cells   !
  !-----------!
  write(f9) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(f9) IN_4 // '<DataArray type='//intp          //  &
                    ' Name="connectivity"'            //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
  write(f9) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Fill up an array with cell offsets and save the header only
  cell_offset = 0
  if(plot_inside) then
    do c1 = c_f, c_l
      cell_offset   = cell_offset + abs(Grid % cells_n_nodes(c1))
      offs_save(c1) = cell_offset
    end do
  else
    do c2 = c_f, c_l
      cell_offset   = cell_offset + Grid % cells_n_nodes(c2)
      offs_save(c2) = cell_offset
    end do
  end if
  call Results % Save_Scalar_Int("offsets", plot_inside,  &
                                  offs_save(c_f:c_l),     &
                                  f8, f9, data_offset, 1)  ! 1 => header only

  ! Fill up an array with cell types and save the header only
  if(plot_inside) then
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .eq. 8) type_save(c1) = VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c1) .eq. 6) type_save(c1) = VTK_WEDGE
      if(Grid % cells_n_nodes(c1) .eq. 4) type_save(c1) = VTK_TETRA
      if(Grid % cells_n_nodes(c1) .eq. 5) type_save(c1) = VTK_PYRAMID
      if(Grid % cells_n_nodes(c1) .lt. 0) type_save(c1) = VTK_POLYHEDRON
    end do
  else
    do c2 = c_f, c_l
      if(Grid % cells_n_nodes(c2) .eq. 4) then
        type_save(c2) = VTK_QUAD
      else if(Grid % cells_n_nodes(c2) .eq. 3) then
        type_save(c2) = VTK_TRIANGLE
      else
        type_save(c2) = VTK_POLYGON
      end if
    end do
  end if
  call Results % Save_Scalar_Int("types", plot_inside,  &
                                  type_save(c_f:c_l),   &
                                  f8, f9, data_offset, 1)  ! 1 => header only

  ! Write parts of header for polyhedral cells
  if(n_polyg > 0) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next


    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faceoffsets"'            //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + (c_l-c_f+1) * IP  ! prepare for next

  end if

  write(f9) IN_3 // '</Cells>'     // LF

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_3 // '<PCellData>' // LF
  end if
  write(f9) IN_3 // '<CellData>' // LF

!------------------------------------------------------------------------------!

  !----------------!
  !                !
  !   Two sweeps   !
  !                !
  !----------------!
  do run = 1, 2

    !------------------------------------------!
    !   Save remnants of the Grid definition   !
    !------------------------------------------!
    if(run .eq. 2) then

      ! Save the nodes' coordinates
      data_size = int(Grid % n_nodes * RP * 3, SP)
      write(f9) data_size
      do n = 1, Grid % n_nodes
        write(f9) Grid % xn(n), Grid % yn(n), Grid % zn(n)
      end do

      ! Save connections
      data_size = int(n_conns * IP, SP)
      write(f9) data_size
      if(plot_inside) then
        do c1 = c_f, c_l

          ! Tetrahedral, pyramid, wedge and hexahedral cells
          if( any( Grid % cells_n_nodes(c1) .eq. (/4,5,6,8/)  ) ) then
            write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c1), c1))-1

          ! Polyhedral cells
          else if(Grid % cells_n_nodes(c1) .lt. 0) then
            write(f9) (Grid % cells_n(1:-Grid % cells_n_nodes(c1), c1))-1

          end if
        end do
      else  ! plot only boundary
        do c2 = c_f, c_l

          ! All cell types
          write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c2), c2))-1

        end do
      end if

      ! Save cell offsets
      call Results % Save_Scalar_Int("offsets", plot_inside,     &
                                      offs_save(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      ! Save cell types
      call Results % Save_Scalar_Int("types", plot_inside,       &
                                      type_save(c_f:c_l),        &
                                      f8, f9, data_offset, run)

      if(n_polyg > 0) then

        ! Write polyhedral cells' faces
        data_size = int(n_polyg * IP, SP)
        write(f9) data_size

        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            write(f9) Grid % cells_n_faces(c1)      ! write number of polygons
            do i_fac = 1, Grid % cells_n_faces(c1)  ! write nodes of each polygn
              s = Grid % cells_f(i_fac, c1)
              if(Grid % faces_s(s) .ne. 0) then  ! face has a shadow, if it ...
                s1 = s                           ! ... is closer, plot that!
                s2 = Grid % faces_s(s)
                dist1 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s1), Grid % yf(s1), Grid % zf(s1))
                dist2 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
                if(dist1 < dist2) s = s1
                if(dist2 < dist1) s = s2
              end if
              n = Grid % faces_n_nodes(s)
              write(f9) n, (Grid % faces_n(1:n, s))-1
            end do
          end if
        end do

        ! Write polyhedral cells' faces offsets
        data_size = int((c_l-c_f+1) * IP, SP)
        write(f9) data_size

        cell_offset = 0
        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            cell_offset = cell_offset + 1           ! to store number of polygs
            do i_fac = 1, Grid % cells_n_faces(c1)  ! to store all the nodes
              s = Grid % cells_f(i_fac, c1)         ! of each polygon
              n = Grid % faces_n_nodes(s)
              cell_offset = cell_offset + 1 + n
            end do
            write(f9) cell_offset

          ! Not a polyhedron, offsets are not needed and set to -1
          else
            write(f9) -1
          end if

        end do

      end if  ! n_polyg > 0

    end if

    !--------------------!
    !   Processor i.d.   !
    !--------------------!
    do c1 = c_f, c_l
      int_save(c1) = Grid % Comm % cell_proc(c1)
    end do
    do c2 = c_f, c_l
      int_save(c2) = Grid % Comm % cell_proc(c2)
    end do
    call Results % Save_Scalar_Int("Grid Processor [1]", plot_inside,   &
                                    int_save(c_f:c_l),                  &
                                    f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(c_f:c_l) = domain
      call Results % Save_Scalar_Int("Grid Domain [1]", plot_inside,  &
                                      int_save(c_f:c_l),              &
                                      f8, f9, data_offset, run)
    end if

!------------------------------------------------------------------------------!
!                         Snapshots area by Yves (2/2)                         !
!------------------------------------------------------------------------------!

! Get the subgrid viscosity when Smagorinsky or Dynamic Smagorinsky or WALE models are used
    kin_vis_t(:) = 0.0
    if(Turb % model .eq. LES_SMAGORINSKY .or.  &
       Turb % model .eq. LES_DYNAMIC .or. &
       Turb % model .eq. LES_WALE) then
       kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l)       
    end if
      
! Sub-Snapshots creation starts here (one for each processor)   
  if(run .eq. 1 .and. plot_inside) then            					     
  			                                                                              													
    if (this_proc < 10)                       write (filename, "(A8,I1,A4)") "SnapShot", this_proc, '.txt'		     
    if (this_proc > 9 .and. this_proc < 100)  write (filename, "(A8,I2,A4)") "SnapShot", this_proc, '.txt'
    if (this_proc > 99)                       write (filename, "(A8,I3,A4)") "SnapShot", this_proc, '.txt' 
  
    open(unit=606+this_proc,file = trim(filename),form='formatted',status='unknown')		     
  
    do i = c_f, c_l									             
       
       if (Grid % vol(i) > 0.000000000001) then							     
          
    	     if(Turb % model .eq. LES_SMAGORINSKY .or.  &
       		Turb % model .eq. LES_DYNAMIC .or. &
		Turb % model .eq. LES_WALE) then
             
        	write (606+this_proc, *) Grid % xc(i), &		 			               
                                 	 Grid % yc(i), &		                                       
        				 Grid % zc(i), &		                                       
        			 	 Flow % u % n(i), &		                                    
        			 	 Flow % v % n(i), &		                                      
        			 	 Flow % w % n(i), &		                                     
        			 	 Flow % p % n(i), &
        			 	 Flow % u % x(i), &		     
        				 Flow % u % y(i), &		      
        				 Flow % u % z(i), &		    
        				 Flow % v % x(i), &		      
        				 Flow % v % y(i), &		   
        				 Flow % v % z(i), &		      
        				 Flow % w % x(i), &		     
        				 Flow % w % y(i), &		    
        				 Flow % w % z(i), &
        				 kin_vis_t(i)		           			 	 
        			 	 
             else
        			 	 
        	write (606+this_proc, *) Grid % xc(i), &		 			               
                                 	 Grid % yc(i), &		                                       
        				 Grid % zc(i), &		                                       
        			 	 Flow % u % n(i), &		                                    
        			 	 Flow % v % n(i), &		                                      
        			 	 Flow % w % n(i), &		                                     
        			 	 Flow % p % n(i), &
              			 	 Flow % u % x(i), &		     
        				 Flow % u % y(i), &		      
        				 Flow % u % z(i), &		    
        				 Flow % v % x(i), &		      
        				 Flow % v % y(i), &		   
        				 Flow % v % z(i), &		      
        				 Flow % w % x(i), &		     
        				 Flow % w % y(i), &		    
        				 Flow % w % z(i)
        				       			 
             end if
        			 	 		 
       end if										            
    end do										            
    
    close(606+this_proc)								            
    
  end if											     

  ! Creation of the total Snapshot			     
  if(run .eq. 2 .and. this_proc .eq. n_proc .and. .not. plot_inside) then	   	             
												     											
    allocate(Snap % x_snap(n_of_cells));     Snap % x_snap = 0.0					     
    allocate(Snap % y_snap(n_of_cells));     Snap % y_snap = 0.0					     
    allocate(Snap % z_snap(n_of_cells));     Snap % z_snap = 0.0
    allocate(Snap % u_snap(n_of_cells));     Snap % u_snap = 0.0					    
    allocate(Snap % v_snap(n_of_cells));     Snap % v_snap = 0.0					     
    allocate(Snap % w_snap(n_of_cells));     Snap % w_snap = 0.0					     
    allocate(Snap % p_snap(n_of_cells));     Snap % p_snap = 0.0					         
    allocate(Snap % dudx_snap(n_of_cells));  Snap % dudx_snap = 0.0					  
    allocate(Snap % dudy_snap(n_of_cells));  Snap % dudy_snap = 0.0					  
    allocate(Snap % dudz_snap(n_of_cells));  Snap % dudz_snap = 0.0					 
    allocate(Snap % dvdx_snap(n_of_cells));  Snap % dvdx_snap = 0.0					 
    allocate(Snap % dvdy_snap(n_of_cells));  Snap % dvdy_snap = 0.0					 
    allocate(Snap % dvdz_snap(n_of_cells));  Snap % dvdz_snap = 0.0					  
    allocate(Snap % dwdx_snap(n_of_cells));  Snap % dwdx_snap = 0.0					  
    allocate(Snap % dwdy_snap(n_of_cells));  Snap % dwdy_snap = 0.0					 
    allocate(Snap % dwdz_snap(n_of_cells));  Snap % dwdz_snap = 0.0
    allocate(Snap % ni_turb(n_of_cells));    Snap % ni_turb = 0.0   					     
  
    do j = 1, n_proc									             
    
      if (j < 10)               write (filename, "(A8,I1,A4)") "SnapShot", j, '.txt'				     
      if (j > 9 .and. j < 100)  write (filename, "(A8,I2,A4)") "SnapShot", j, '.txt'
      if (j > 99)               write (filename, "(A8,I3,A4)") "SnapShot", j, '.txt'
  
      open(unit=606+j,file = trim(filename),form='formatted',status='unknown')			     

      do
      
    	if(Turb % model .eq. LES_SMAGORINSKY .or.  &
      	   Turb % model .eq. LES_DYNAMIC .or. &
	   Turb % model .eq. LES_WALE) then                                                          
      					
	    	read (606+j,*, end=10)  Snap % x_snap(n_c_tot+1), &				             
				   	Snap % y_snap(n_c_tot+1), &				             
				   	Snap % z_snap(n_c_tot+1), &				             
				   	Snap % u_snap(n_c_tot+1), &				             
				   	Snap % v_snap(n_c_tot+1), &					     
				   	Snap % w_snap(n_c_tot+1), &				             
				   	Snap % p_snap(n_c_tot+1), &
				        Snap % dudx_snap(n_c_tot+1), &					     
				   	Snap % dudy_snap(n_c_tot+1), &					 
				   	Snap % dudz_snap(n_c_tot+1), &					   
				   	Snap % dvdx_snap(n_c_tot+1), &					     
				   	Snap % dvdy_snap(n_c_tot+1), &					      
				   	Snap % dvdz_snap(n_c_tot+1), &					   
				   	Snap % dwdx_snap(n_c_tot+1), &					      
				   	Snap % dwdy_snap(n_c_tot+1), &					   
				   	Snap % dwdz_snap(n_c_tot+1), &				   	
				   	Snap % ni_turb(n_c_tot+1)
				   
	    else
	    
	    	read (606+j,*, end=10)  Snap % x_snap(n_c_tot+1), &				             
				   	Snap % y_snap(n_c_tot+1), &				             
				   	Snap % z_snap(n_c_tot+1), &				             
				   	Snap % u_snap(n_c_tot+1), &				             
				   	Snap % v_snap(n_c_tot+1), &					     
				   	Snap % w_snap(n_c_tot+1), &				             
				   	Snap % p_snap(n_c_tot+1), &
				        Snap % dudx_snap(n_c_tot+1), &					     
				   	Snap % dudy_snap(n_c_tot+1), &					 
				   	Snap % dudz_snap(n_c_tot+1), &					   
				   	Snap % dvdx_snap(n_c_tot+1), &					     
				   	Snap % dvdy_snap(n_c_tot+1), &					      
				   	Snap % dvdz_snap(n_c_tot+1), &					   
				   	Snap % dwdx_snap(n_c_tot+1), &					      
				   	Snap % dwdy_snap(n_c_tot+1), &					   
				   	Snap % dwdz_snap(n_c_tot+1)			   
				   
            end if	   
				   				   				   				       			   
	    n_c_tot = n_c_tot + 1								     
	    
      end do										             

      10 close (606+j, status='delete')								
      								      								      								    
    end do
    
    ! Ordering of the Snapshot according to cells' coordinates											     
    allocate(Snap % pos_ind(n_c_tot));  Snap % pos_ind = 0    
    allocate(Snap % x_2(n_c_tot));      Snap % x_2 = 0.0
    allocate(Snap % y_2(n_c_tot));      Snap % y_2 = 0.0
    allocate(Snap % z_2(n_c_tot));      Snap % z_2 = 0.0

    Snap % x_2(:) = HUGE
    Snap % y_2(:) = HUGE
    
    do i = 1, n_c_tot
    
      Snap % pos_ind(i) = i
      Snap % x_2(i) = min(Snap % x_2(i),  &
                     Snap % x_snap(i))
      Snap % y_2(i) = min(Snap % y_2(i),  &
                     Snap % y_snap(i))
      Snap % z_2(i) = Snap % z_snap(i)
      
    end do
    
  call Sort % Three_Real_Carry_Int(Snap % x_2, Snap % y_2, Snap % z_2, Snap % pos_ind)		          										
    
    ! Creation of the output Snapshot .txt file
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

    open(unit=607,file = trim(filename),form='formatted',status='unknown')				  
    
    do i = 1, n_c_tot
      
    	if(Turb % model .eq. LES_SMAGORINSKY .or.  &
           Turb % model .eq. LES_DYNAMIC .or. &
	   Turb % model .eq. LES_WALE) then											  
      
     		 write (607, *) Snap % x_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % y_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % z_snap(Snap % pos_ind(i)), ',',  &						    
          	     		Snap % u_snap(Snap % pos_ind(i)), ',',  &						   
        	     		Snap % v_snap(Snap % pos_ind(i)), ',',  &						   
        	     		Snap % w_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % p_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudz_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdz_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdz_snap(Snap % pos_ind(i)), ',',  &        	     		
        	     		Snap % ni_turb(Snap % pos_ind(i))
        	     		
            else
            
     		 write (607, *) Snap % x_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % y_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % z_snap(Snap % pos_ind(i)), ',',  &						    
          	     		Snap % u_snap(Snap % pos_ind(i)), ',',  &						   
        	     		Snap % v_snap(Snap % pos_ind(i)), ',',  &						   
        	     		Snap % w_snap(Snap % pos_ind(i)), ',',  &						    
        	     		Snap % p_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dudz_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dvdz_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdx_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdy_snap(Snap % pos_ind(i)), ',',  &
        	     		Snap % dwdz_snap(Snap % pos_ind(i))
            
            end if							           							
        								     
    end do
    
    close(607)												  
    
    deallocate(Snap % u_snap)										
    deallocate(Snap % v_snap)										  
    deallocate(Snap % w_snap)										  
    deallocate(Snap % p_snap)										 
    deallocate(Snap % x_snap)										  
    deallocate(Snap % y_snap)										  
    deallocate(Snap % z_snap)										  
    deallocate(Snap % pos_ind)										  
    deallocate(Snap % x_2)										 
    deallocate(Snap % y_2)										  
    deallocate(Snap % z_2)
    deallocate(Snap % dudx_snap)
    deallocate(Snap % dudy_snap)
    deallocate(Snap % dudz_snap)
    deallocate(Snap % dvdx_snap)
    deallocate(Snap % dvdy_snap)
    deallocate(Snap % dvdz_snap)
    deallocate(Snap % dwdx_snap)
    deallocate(Snap % dwdy_snap)
    deallocate(Snap % dwdz_snap)
    deallocate(Snap % ni_turb)										  
         												 
    n_c_tot = 0;											 

  end if
  
!------------------------------------------------------------------------------!												  
  

    !--------------!
    !   Velocity   !
    !--------------!
    call Results % Save_Vector_Real("Velocity [m/s]", plot_inside,  &
                                    Flow % u % n(c_f:c_l),          &
                                    Flow % v % n(c_f:c_l),          &
                                    Flow % w % n(c_f:c_l),          &
                                    f8, f9, data_offset, run)
    !--------------------!
    !   Courant number   !
    !--------------------!
    if(plot_inside) then
      call Flow % Calculate_Courant_In_Cells(save_01)
      call Results % Save_Scalar_Real("Courant Number [1]", plot_inside,  &
                                      save_01(c_f:c_l),                   &
                                      f8, f9, data_offset, run)
    end if

    !---------------!
    !   Potential   !
    !---------------!
    call Results % Save_Scalar_Real("Potential [m^2/s]", plot_inside,  &
                                    Flow % pot % n(c_f:c_l),           &
                                    f8, f9, data_offset, run)

    !--------------------------------------!
    !   Pressure correction and pressure   !
    !--------------------------------------!
    call Results % Save_Scalar_Real("Pressure Correction [Pa]",  &
                                    plot_inside,                 &
                                    Flow % pp % n(c_f:c_l),      &
                                    f8, f9, data_offset, run)
    save_01(:) = 0.0
    save_02(:) = 0.0
    save_03(:) = 0.0
    do c1 = c_f, c_l
      save_01(c1) = Flow % pp % x(c1) * Grid % vol(c1)
      save_02(c1) = Flow % pp % y(c1) * Grid % vol(c1)
      save_03(c1) = Flow % pp % z(c1) * Grid % vol(c1)
    end do
    call Results % Save_Vector_Real("Pressure Correction Force [N]",  &
                                    plot_inside,                      &
                                    save_01(c_f:c_l),                 &
                                    save_02(c_f:c_l),                 &
                                    save_03(c_f:c_l),                 &
                                    f8, f9, data_offset, run)

    call Results % Save_Scalar_Real("Pressure [Pa]", plot_inside,  &
                                    Flow % p % n(c_f:c_l),         &
                                    f8, f9, data_offset, run)
    save_01(:) = 0.0
    save_02(:) = 0.0
    save_03(:) = 0.0
    do c1 = c_f, c_l
      save_01(c1) = Flow % p % x(c1) * Grid % vol(c1)
      save_02(c1) = Flow % p % y(c1) * Grid % vol(c1)
      save_03(c1) = Flow % p % z(c1) * Grid % vol(c1)
    end do
    call Results % Save_Vector_Real("PressureForce [N]", plot_inside,    &
                                    save_01(c_f:c_l),                    &
                                    save_02(c_f:c_l),                    &
                                    save_03(c_f:c_l),                    &
                                    f8, f9, data_offset, run)

    !-----------------!
    !   Temperature   !
    !-----------------!
    if(Flow % heat_transfer) then
      call Results % Save_Scalar_Real("Temperature [K]", plot_inside,  &
                                      Flow % t % n(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0
      do c1 = c_f, c_l
        save_01(c1) = Flow % t % x(c1)
        save_02(c1) = Flow % t % y(c1)
        save_03(c1) = Flow % t % z(c1)
      end do

      if(.not. Flow % mass_transfer) then
        call Flow % Grad_Variable(Flow % t)
      else
        call Vof % Calculate_Grad_Matrix_With_Front()
        call Vof % Grad_Variable_With_Front(Flow % t, Vof % t_sat)
        call Flow % Calculate_Grad_Matrix()
      end if

      call Results % Save_Vector_Real("Temperature Gradients [K/m]",  &
                                      plot_inside,                    &
                                      save_01(c_f:c_l),               &
                                      save_02(c_f:c_l),               &
                                      save_03(c_f:c_l),               &
                                      f8, f9, data_offset, run)

    end if

    !-------------------------!
    !   Physical properties   !
    !-------------------------!
    call Results % Save_Scalar_Real("Physical Density [kg/m^3]",      &
                                    plot_inside,                      &
                                    Flow % density(c_f:c_l),          &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Physical Viscosity [Pa s]",      &
                                    plot_inside,                      &
                                    Flow % viscosity(c_f:c_l),        &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Physical Conductivity [W/m/K]",  &
                                    plot_inside,                      &
                                    Flow % conductivity(c_f:c_l),     &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Physical Capacity [J/K]",        &
                                    plot_inside,                      &
                                    Flow % capacity(c_f:c_l),         &
                                    f8, f9, data_offset, run)

    if(Turb % rough_walls) then
      call Results % Save_Scalar_Real("Roughness Coefficient z_o [1]",  &
                                      plot_inside,                      &
                                      Turb % z_o(c_f:c_l),              &
                                      f8, f9, data_offset, run)

    end if

    !---------------------!
    !   Volume fraction   !
    !---------------------!
    if(Flow % with_interface) then
      call Results % Save_Scalar_Real("Vof Sharp [1]",                  &
                                      plot_inside,                      &
                                      Vof % fun % n(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Vof Smooth [1]",                 &
                                      plot_inside,                      &
                                      Vof % smooth % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Vof Curvature [1/m]",            &
                                      plot_inside,                      &
                                      Vof % curv(c_f:c_l),              &
                                      f8, f9, data_offset, run)
      call Results % Save_Vector_Real("Vof SurfaceNormals [1]",         &
                                      plot_inside,                      &
                                      Vof % nx(c_f:c_l),                &
                                      Vof % ny(c_f:c_l),                &
                                      Vof % nz(c_f:c_l),                &
                                      f8, f9, data_offset, run)
      call Results % Save_Vector_Real("Vof SurfaceTensionForce [N]",    &
                                      plot_inside,                      &
                                      Vof % surf_fx(c_f:c_l),           &
                                      Vof % surf_fy(c_f:c_l),           &
                                      Vof % surf_fz(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      if (allocated(Vof % m_dot)) then
        call Results % Save_Scalar_Real("Vof MassTransfer [kg/m^3/s]",  &
                                        plot_inside,                    &
                                        Vof % m_dot(c_f:c_l),           &
                                        f8, f9, data_offset, run)
      end if
    end if

    !---------------------------------------!
    !   Number of impacts and reflections   !
    !---------------------------------------!
    if(Flow % with_particles .and. .not. plot_inside) then
      call Results % Save_Scalar_Real("Particles Reflected [1]",     &
                                      plot_inside,                   &
                                      Swarm % n_reflected(c_f:c_l),  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Particles Deposited [1]",     &
                                      plot_inside,                   &
                                      Swarm % n_deposited(c_f:c_l),  &
                                      f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      call Results % Save_Scalar_Real(phi % name, plot_inside,   &
                                      phi % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
    end do

    !-----------------!
    !   Q-criterion   !
    !-----------------!
    call Calculate_Shear_And_Vorticity(Flow)
    phi_save(:) = 0.0
    do c1 = c_f, c_l
      phi_save(c1) = (Flow % vort(c1)**2 - Flow % shear(c1)**2)/4.
    end do
    call Results % Save_Scalar_Real("Q Criterion [1/s^2]", plot_inside,   &
                                    phi_save(c_f:c_l),                    &
                                    f8, f9, data_offset, run)

    !--------------------------!
    !   Turbulent quantities   !
    !--------------------------!

    ! Save kin and eps
    if(Turb % model .eq. K_EPS                 .or.  &
       Turb % model .eq. K_EPS_ZETA_F          .or.  &
       Turb % model .eq. HYBRID_LES_RANS       .or.  &
       Turb % model .eq. RSM_MANCEAU_HANJALIC  .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC  ) then
      call Results % Save_Scalar_Real("Turbulent Kinetic Energy [m^2/s^2]",  &
                            plot_inside,                                     &
                            Turb % kin % n(c_f:c_l),                         &
                            f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Turbulent Dissipation [m^2/s^3]",    &
                                      plot_inside,                          &
                                      Turb % eps % n(c_f:c_l),              &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real(                                        &
                            "Turbulent Kinetic Energy Production [m^2/s^3]",  &
                            plot_inside,                                      &
                            Turb % p_kin(c_f:c_l),                            &
                            f8, f9, data_offset, run)
    end if

    ! Save zeta and f22
    if(Turb % model .eq. K_EPS_ZETA_F .or.  &
       Turb % model .eq. HYBRID_LES_RANS) then
      v2_calc(:) = 0.0
      do c1 = c_f, c_l
        v2_calc(c1) = Turb % kin % n(c1) * Turb % zeta % n(c1)
      end do
!      call Results % Save_Scalar_Real("Turbulent Quantity V2 [m^2/s^2]",    &
!                                      plot_inside,                          &
!                                      v2_calc (c_f:c_l),                    &
!                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Turbulent Quantity Zeta [1]",        &
                                      plot_inside,                          &
                                      Turb % zeta % n(c_f:c_l),             &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Turbulent Quantity F22 [1]",         &
                                      plot_inside,                          &
                                      Turb % f22  % n(c_f:c_l),             &
                                      f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("Turbulent Quantity T2 [K^2]",      &
                                        plot_inside,                        &
                                        Turb % t2 % n(c_f:c_l),             &
                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("Turbulent T2 Production [K^2/s]",  &
!                                        plot_inside,                        &
!                                        Turb % p_t2(c_f:c_l),               &
!                                        f8, f9, data_offset, run)
        call Results % Save_Vector_Real("Turbulent Heat Flux [K m/s]",  &
                                        plot_inside,                    &
                                        Turb % ut % n(c_f:c_l),         &
                                        Turb % vt % n(c_f:c_l),         &
                                        Turb % wt % n(c_f:c_l),         &
                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("Turbulent Quantity Alpha L",     &
!                                        plot_inside,                      &
!                                        Turb % alpha_l(c_f:c_l),          &
!                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("Turbulent Quantity Alpha U",     &
!                                        plot_inside,                      &
!                                        Turb % alpha_u(c_f:c_l),          &
!                                        f8, f9, data_offset, run)
      end if
    end if

    if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Results % Save_Scalar_Real("Turbulent Quantity F22 [1]",  &
                                      plot_inside,                   &
                                      Turb % f22 % n(c_f:c_l),       &
                                      f8, f9, data_offset, run)
    end if

    ! Save vis and vis_t
    if(Turb % model .eq. DES_SPALART .or.  &
       Turb % model .eq. SPALART_ALLMARAS) then
      call Results % Save_Scalar_Real("Turbulent Viscosity [Pa s]",  &
                                      plot_inside,                   &
                                      Turb % vis % n(c_f:c_l),       &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Vorticity Magnitude [1/s]",   &
                                      plot_inside,                   &
                                      Flow % vort(c_f:c_l),          &
                                      f8, f9, data_offset, run)
    end if

    kin_vis_t(:) = 0.0
    if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       Turb % model .ne. HYBRID_LES_RANS     .and.  &
       Turb % model .ne. DNS) then
      kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real("Eddy Over Molecular Viscosity [1]",  &
                                      plot_inside,                          &
                                      kin_vis_t(c_f:c_l),                   &
                                      f8, f9, data_offset, run)
    end if

    if(Turb % model .eq. HYBRID_LES_RANS) then
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real(                                       &
                                  "Rans Eddy Over Molecular Viscosity [1]",  &
                                  plot_inside,                               &
                                  kin_vis_t(c_f:c_l),                        &
                                  f8, f9, data_offset, run)
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = Turb % vis_t_sgs(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real(                                      &
                                  "Sgs Eddy Over Molecular Viscosity [1]",  &
                                  plot_inside,                              &
                                  kin_vis_t(c_f:c_l),                       &
                                  f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

      ! Note: follows the order in which Paraview stores tensors
      call Results % Save_Tensor_6_Real("Reynolds Stress [m^2/s^2]",  &
                                        plot_inside,                  &
                                        Turb % uu % n(c_f:c_l),       &
                                        Turb % vv % n(c_f:c_l),       &
                                        Turb % ww % n(c_f:c_l),       &
                                        Turb % uv % n(c_f:c_l),       &
                                        Turb % vw % n(c_f:c_l),       &
                                        Turb % uw % n(c_f:c_l),       &
                                        f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Vector_Real("Turbulent Heat Flux [K m/s]",  &
                                        plot_inside,                    &
                                        Turb % ut % n(c_f:c_l),         &
                                        Turb % vt % n(c_f:c_l),         &
                                        Turb % wt % n(c_f:c_l),         &
                                        f8, f9, data_offset, run)
      end if
    end if

    ! Statistics for large-scale simulations of turbulence
    if(Turb % statistics) then
      call Results % Save_Vector_Real("Mean Velocity [m/s]",      &
                                      plot_inside,                &
                                      Turb % u_mean(c_f:c_l),     &
                                      Turb % v_mean(c_f:c_l),     &
                                      Turb % w_mean(c_f:c_l),     &
                                      f8, f9, data_offset, run)
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0
      save_04(:) = 0.0
      save_05(:) = 0.0
      save_06(:) = 0.0

      ! Note: follows the order in which Paraview stores tensors
      do c1 = c_f, c_l
        save_01(c1) = Turb % uu_res(c1) - Turb % u_mean(c1) * Turb % u_mean(c1)
        save_02(c1) = Turb % vv_res(c1) - Turb % v_mean(c1) * Turb % v_mean(c1)
        save_03(c1) = Turb % ww_res(c1) - Turb % w_mean(c1) * Turb % w_mean(c1)
        save_04(c1) = Turb % uv_res(c1) - Turb % u_mean(c1) * Turb % v_mean(c1)
        save_05(c1) = Turb % vw_res(c1) - Turb % v_mean(c1) * Turb % w_mean(c1)
        save_06(c1) = Turb % uw_res(c1) - Turb % u_mean(c1) * Turb % w_mean(c1)
      end do
      call Results % Save_Tensor_6_Real("Mean Reynolds Stress [m^s/s^2]",  &
                                        plot_inside,                       &
                                        save_01(c_f:c_l),                  &
                                        save_02(c_f:c_l),                  &
                                        save_03(c_f:c_l),                  &
                                        save_04(c_f:c_l),                  &
                                        save_05(c_f:c_l),                  &
                                        save_06(c_f:c_l),                  &
                                        f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("Mean Temperature [K]",           &
                                        plot_inside,                      &
                                        Turb % t_mean(c_f:c_l),           &
                                        f8, f9, data_offset, run)
        phi_save(:) = 0.0
        save_01(:) = 0.0
        save_02(:) = 0.0
        save_03(:) = 0.0
        do c1 = c_f, c_l
          phi_save(c1) = Turb % t2_res(c1) - Turb % t_mean(c1)*Turb % t_mean(c1)
          save_01(c1)  = Turb % ut_res(c1) - Turb % u_mean(c1)*Turb % t_mean(c1)
          save_02(c1)  = Turb % vt_res(c1) - Turb % v_mean(c1)*Turb % t_mean(c1)
          save_03(c1)  = Turb % wt_res(c1) - Turb % w_mean(c1)*Turb % t_mean(c1)
        end do
        call Results % Save_Scalar_Real("Mean Turbulent Quantity T2 [K^2]",    &
                                        plot_inside,                           &
                                        phi_save(c_f:c_l),                     &
                                        f8, f9, data_offset, run)
        call Results % Save_Vector_Real("Mean Turbulent Heat Flux [K m/s]",  &
                                        plot_inside,                         &
                                        save_01(c_f:c_l),                    &
                                        save_02(c_f:c_l),                    &
                                        save_03(c_f:c_l),                    &
                                        f8, f9, data_offset, run)
      end if

      ! Scalars
      do sc = 1, Flow % n_scalars
        phi => Flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        do c1 = c_f, c_l
          phi_save(c1) = Turb % scalar_mean(sc, c1)
        end do
        call Results % Save_Scalar_Real(name_mean, plot_inside,    &
                                        phi_save(c_f:c_l),         &
                                        f8, f9, data_offset, run)
      end do
    end if

    ! Save y+ for all turbulence models and ILES (none option)
    if(Turb % model .ne. DNS) then
      call Results % Save_Scalar_Real("Turbulent Quantity Y Plus [1]",  &
                                      plot_inside,                      &
                                      Turb % y_plus(c_f:c_l),           &
                                      f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    call Results % Save_Scalar_Real("Grid Cell Volume [m^3]",      &
                                    plot_inside,                   &
                                    Grid % vol(c_f:c_l),           &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Grid Wall Distance [m]",      &
                                    plot_inside,                   &
                                    Grid % wall_dist(c_f:c_l),     &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Grid Cell Delta Max [m]",     &
                                    plot_inside,                   &
                                    Turb % h_max(c_f:c_l),         &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Grid Cell Delta Min [m]",     &
                                    plot_inside,                   &
                                    Turb % h_min(c_f:c_l),         &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("Grid Cell Delta Wall [m]",    &
                                    plot_inside,                   &
                                    Turb % h_w  (c_f:c_l),         &
                                    f8, f9, data_offset, run)

    !---------------------------------------------------------------------!
    !   Variables in the first computational point, plotted at boundary   !
    !---------------------------------------------------------------------!

    ! Engage only for boundary plots (not inside means on the boundary)
    if( .not. plot_inside ) then 

      ! Initialize working variables to zero
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0

      ! Copy internal values to boundary
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then 
          save_01(c2) = Flow % u % n(c1)
          save_02(c2) = Flow % v % n(c1)
          save_03(c2) = Flow % w % n(c1)
        end if
      end do
      call Results % Save_Vector_Real("Velocity Near Wall [m/s]", plot_inside, &
                                      save_01(c_f:c_l),                        &
                                      save_02(c_f:c_l),                        &
                                      save_03(c_f:c_l),                        &
                                      f8, f9, data_offset, run)

      if(Turb % model .eq. K_EPS                 .or.  &
         Turb % model .eq. K_EPS_ZETA_F          .or.  &
         Turb % model .eq. HYBRID_LES_RANS) then

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do s = 1, Grid % n_faces
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          if(c2 < 0) then
            var_ins(c2) = Turb % kin % n(c1)
          end if
        end do

        call Results % Save_Scalar_Real("T.K.E. Near Wall [m^2/s^2]",  &
                                        plot_inside,                   &
                                        var_ins(c_f:c_l),              &
                                        f8, f9, data_offset, run)

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do s = 1, Grid % n_faces
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          if(c2 < 0) then
            var_ins(c2) = Turb % y_plus(c1)
          end if
        end do

        call Results % Save_Scalar_Real("y+ Near Wall [1]",         &
                                        plot_inside,                &
                                        var_ins(c_f:c_l),           &
                                        f8, f9, data_offset, run)

        do sc = 1, Flow % n_scalars
          phi => Flow % scalar(sc)
          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi % n(c1)
            end if
          end do

          call Results % Save_Scalar_Real("Scalar Near Wall",        &
                                          plot_inside,               &
                                          var_ins(c_f:c_l),          &
                                          f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi_save(c1)  ! Turb % scalar_mean(sc, c1)
            end if
          end do

          call Results % Save_Scalar_Real("Mean Scalar Near Wall",  &
                                plot_inside,                        &
                                var_ins(c_f:c_l),                   &
                                f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi % q(c2)  ! Turb % scalar_mean(sc, c1)
            end if
          end do

          call Results % Save_Scalar_Real("Wall Scalar Flux",        &
                                          plot_inside,               &
                                          var_ins(c_f:c_l),          &
                                          f8, f9, data_offset, run)

        end do
      end if
    end if

    !----------------------!
    !   Save user arrays   !
    !----------------------!
    do ua = 1, Grid % n_user_arrays

      a_name = 'A_00'
      write(a_name(3:4), '(i2.2)') ua
      call Results % Save_Scalar_Real(a_name, plot_inside,            &
                                      Grid % user_array(ua,c_f:c_l),  &
                                      f8, f9, data_offset, run)
    end do

    !----------------------!
    !                      !
    !   End of cell data   !
    !                      !
    !----------------------!
    if(run .eq. 1) then
      if(n_proc > 1 .and. this_proc .eq. 1) then
        write(f8) IN_3 // '</PCellData>' // LF
      end if
      write(f9) IN_3 // '</CellData>' // LF

      write(f9) IN_2 // '</Piece>'            // LF
      write(f9) IN_1 // '</UnstructuredGrid>' // LF

      !-------------------!
      !                   !
      !   Appended data   !
      !                   !
      !-------------------!
      write(f9) IN_0 // '<AppendedData encoding="raw">' // LF
      write(f9) '_'
    else
      write(f9) LF // IN_0 // '</AppendedData>' // LF
    end if

  end do  ! run
  
!------------------------------------------------------------------------------!

  write(f9) IN_0 // '</VTKFile>'          // LF
  close(f9)

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1) then
    do n = 1, n_proc
      if(plot_inside) then
        call File % Set_Name(name_out_9,        &
                             processor=n,       &
                             time_step=ts,      &
                             extension='.vtu',  &
                             domain=domain)
      else
        call File % Set_Name(name_out_9,        &
                             processor=n,       &
                             time_step=ts,      &
                             appendix ='-bnd',  &
                             extension='.vtu',  &
                             domain=domain)
      end if
      write(f8) IN_2 // '<Piece Source="', trim(name_out_9), '"/>' // LF
    end do
    write(f8) IN_1 // '</PUnstructuredGrid>' // LF
    write(f8) IN_0 // '</VTKFile>'           // LF
    close(f8)
  end if

  call Work % Disconnect_Int_Cell(int_save, type_save, offs_save)
  call Work % Disconnect_Real_Cell(save_01, save_02, save_03,  &
                                   save_04, save_05, save_06)
  call Work % Disconnect_Real_Cell(var_ins, v2_calc, kin_vis_t, phi_save)

  call Profiler % Stop('Save_Vtu_Results')

  end subroutine
