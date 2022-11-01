!==============================================================================!
  subroutine Save_Results(Results,  &
                          Flow, Turb, Vof, Swarm, ts, plot_inside, domain)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
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
  
  
  
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------DANGER: STARTING YVES ZONE (1/2)
  
  type SnapShot	! Hay que definir un nuevo Type, con vectores para las coordenadas, la velocidad, y la presión. Cuando tengamos recogida la información para todos los 	
  		! subdominios (en Files distintos, uno por cada procesador), UNO de los procesadores se encargará de recopilarla toda en estos vectores; así, tendremos
  		! la coordenada "x" de CADA celda del dominio en un sólo vector, la "y" en otro, la presión en otro, &c.

    ! Los vectores individuales que queremos:
    real, allocatable :: u_snap(:), v_snap(:), w_snap(:), p_snap(:)		  ! Streamwise, Vertical, and Spanwise VELOCITY in the Cell ; Pressure in the Cell
    real, allocatable :: x_snap(:), y_snap(:), z_snap(:)			  ! "x", "y", and "z" position of the Cell
    
    real, allocatable :: dudx_snap(:), dudy_snap(:), dudz_snap(:)		  ! "x", "y", and "z" GRADIENTS of Streamwise Velocity
    real, allocatable :: dvdx_snap(:), dvdy_snap(:), dvdz_snap(:)		  ! "x", "y", and "z" GRADIENTS of Vertical Velocity
    real, allocatable :: dwdx_snap(:), dwdy_snap(:), dwdz_snap(:)		  ! "x", "y", and "z" GRADIENTS of Spanwise Velocity
    
    integer, allocatable :: puesto(:)			  ! Position of the Cell in the ordered Single Vector
    
    real,    allocatable :: equis2(:)			  ! "x" coordinate of the Cell. ¡¡NECESARY!! Otherwise "z_snap" would be ordered differently after Sorting
    real,    allocatable :: ygriega2(:)			  ! "y" coordinate of the Cell. ¡¡NECESARY!! Otherwise "z_snap" would be ordered differently after Sorting
    real,    allocatable :: ceta2(:)			  ! "z" coordinate of the Cell. ¡¡NECESARY!! Otherwise "z_snap" would be ordered differently after Sorting
    
  end type
  
  type(SnapShot) :: Snap		! Creamos nuestra variable "Snap", del tipo "SnapShot
  
  ! Ahora algunas variables locales:
  integer :: turb_interval2 = 1430	! Intervalo de recogida de SnapShots (INPUT)
  integer :: IniTurb = 253110		! Primer TS en el que se recogen SnapShots (INPUT)
  integer :: ene = 0			! Número de TS que han pasado desde el primer TS en el que se recogen SnapShots
  
  integer :: i = 0			! Un contador
  integer :: j = 0			! Un contador
  
  integer :: n_of_cells = 1200000	! Número de Cells del dominio completo; IMPORTANTE aproximar POR LO ALTO, ya que se usa para alocar espacio en los vectores. (INPUT)
  
  character(len=1024) :: filename	! Nombre del File que queremos abrir o cerrar.
  
  integer :: n_c_tot = 0  		! Contador del número total de Cells; lo primero que hace el procesador ÚNICO es abrir los Files individuales de cada subdominio
  					! e ir construyendo los vectores ÚNICOS con los valores que lée; llevar la cuenta de cuántas Cells llevas leídas es necesario
  					! para poner cada Cell detrás del anterior, especialmente cuando cambias de la última de un subdominio a la primera de otro
  
  integer :: Nceta = 37  		! Número de celdas a lo largo de la dirección homogénea "z"
  					
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ENDING YVES ZONE (1/2)



!------------------------------[Local parameters]------------------------------!
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

  call Cpu_Timer % Start('Save_Vtu_Results')

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



!---------------------------------------------------------------------------------------------------------------------------------------------------------------------DANGER: STARTING YVES ZONE (2/2)

  ene = ts - IniTurb											! Check how many TSs have past since the one where we started to gather SnapShots
  
!----------------------------------------------------------------------------------------------------------------------------------------For SPECTRA

!adsf

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------DANGER: STARTING YVES ZONE (2/2)
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------GENERAL SNAPSHOTS:
  
  if(ene > -1 .and. mod(ene, turb_interval2) < 1 .and. run .eq. 1 .and. plot_inside) then		! If we HAVE started to gather SnapShots, and this is a TS where we're supposed to
  													! save one, do so, but only once (in Run 1) and only for 3D Cells (plot_inside):
  
    if (this_proc < 10) write (filename, "(A8,I1,A4)") "SnapShot", this_proc, '.txt'			  ! The name of the SUBsnapshot file...
    if (this_proc > 9) write (filename, "(A8,I2,A4)") "SnapShot", this_proc, '.txt'			  ! ...is "SnapShot(ProcessorNumber).txt"
  
    open(unit=606+this_proc,file = trim(filename),form='formatted',status='unknown')			  ! So go ahead and create the file...
  
    do i = c_f, c_l											  ! ...and for every Cell of the SubDomain...
      if (Grid % vol(i) > 0.000000000001) then								    ! ...(every 3D Cell, that is)...
        write (606+this_proc, *) Grid % xc(i),  &		 					      ! ...write a row of the file containing its "x" Coordinate,...
        								     Grid % yc(i),  &		      ! ..."y" Coordinate,...
        								     Grid % zc(i),  &		      ! ..."z" Coordinate,...
        								  Flow % u % n(i),  &		      ! ...Streamwise Velocity,...
        								  Flow % v % n(i),  &		      ! ...Vertical Velocity,...
        								  Flow % w % n(i),  &		      ! ...Spanwise Velocity,...
        								  Flow % p % n(i),  &		      ! ...and Pressure...
        								  Flow % u % x(i),  &		      ! asdf
        								  Flow % u % y(i),  &		      ! asdf
        								  Flow % u % z(i),  &		      ! asdf
        								  Flow % v % x(i),  &		      ! asdf
        								  Flow % v % y(i),  &		      ! asdf
        								  Flow % v % z(i),  &		      ! asdf
        								  Flow % w % x(i),  &		      ! asdf
        								  Flow % w % y(i),  &		      ! asdf
        								  Flow % w % z(i)		      ! asdf
      end if												    ! (END "if 3D Cell")
    end do												  ! (END "for every Cell of the SubDomain")
    
    close(606+this_proc)										  ! ...then close the file and go on with the module
    
  end if												! [END "If we HAVE started to gather SnapShots, and this is a TS where (...)]"

  
  
  if(ene > -1 .and. mod(ene, turb_interval2) < 1 .and. run .eq. 2  &					! Now, if we HAVE started to gather SnapShots, and this is a TS where we're supposed to
  	   .and. this_proc .eq. n_proc .and. .not. plot_inside) then					! save one, AND we've already gathered the data about every Cell of every SubDomain
													! ([last part (.not. plot_inside) of the last Run (Run 2)], then use ONLY ONE PROCESSOR
													! (n_proc) to read through the individual files and join all info in a SINGLE SnapShot:
												
    allocate(Snap % u_snap(n_of_cells));  Snap % u_snap = 0.						  !  Start by allocating space for every "Single Vector". We haven't counted...
    allocate(Snap % v_snap(n_of_cells));  Snap % v_snap = 0.						  !  ...exactly how many Cells are in the whole domain yet (because we'll do...
    allocate(Snap % w_snap(n_of_cells));  Snap % w_snap = 0.						  !  ...it WHILE we build the vectors), so for the allocation we'll use the...
    allocate(Snap % p_snap(n_of_cells));  Snap % p_snap = 0.						  !  ...rough estimate that we provided earlier as INPUT; "n_of_cells".
    allocate(Snap % x_snap(n_of_cells));  Snap % x_snap = 0.						  !  ------------------------------------------------------------------
    allocate(Snap % y_snap(n_of_cells));  Snap % y_snap = 0.						  !  ------------------------------------------------------------------
    allocate(Snap % z_snap(n_of_cells));  Snap % z_snap = 0.						  !  ------------------------------------------------------------------
    allocate(Snap % dudx_snap(n_of_cells));  Snap % dudx_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dudy_snap(n_of_cells));  Snap % dudy_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dudz_snap(n_of_cells));  Snap % dudz_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dvdx_snap(n_of_cells));  Snap % dvdx_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dvdy_snap(n_of_cells));  Snap % dvdy_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dvdz_snap(n_of_cells));  Snap % dvdz_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dwdx_snap(n_of_cells));  Snap % dwdx_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dwdy_snap(n_of_cells));  Snap % dwdy_snap = 0.					  !  ------------------------------------------------------------------
    allocate(Snap % dwdz_snap(n_of_cells));  Snap % dwdz_snap = 0.					  !  --------------------------------------------------------------------------
  
    do j = 1, n_proc											  ! Then, for every SUBsnapshot...
    
      if (j < 10) write (filename, "(A8,I1,A4)") "SnapShot", j, '.txt'					    ! ...[whose names, again, are...
      if (j > 9) write (filename, "(A8,I2,A4)") "SnapShot", j, '.txt'					    ! ..."SnapShot(ProcessorNumber).txt"]...
  
      open(unit=606+j,file = trim(filename),form='formatted',status='unknown')				    ! ...go ahead and open the SUBfile:

      do												    ! For every ROW* of the SUBfile, we're going to...
      					
	    read (606+j,*, end=10) Snap % x_snap(n_c_tot+1),     &					      ! ...read the "x" coordinate,...
				   Snap % y_snap(n_c_tot+1),     &					      ! ..."y" coordinate,...
				   Snap % z_snap(n_c_tot+1),     &					      ! ..."z" Coordinate,...
				   Snap % u_snap(n_c_tot+1),     &					      ! ...Streamwise Velocity,...
				   Snap % v_snap(n_c_tot+1),     &					      ! ...Vertical Velocity,...
				   Snap % w_snap(n_c_tot+1),     &					      ! ...Spanwise Velocity,...
				   Snap % p_snap(n_c_tot+1),     &					      ! ...and Pressure, of the corresponding CELL.
				   Snap % dudx_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dudy_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dudz_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dvdx_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dvdy_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dvdz_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dwdx_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dwdy_snap(n_c_tot+1),  &					      ! asdf
				   Snap % dwdz_snap(n_c_tot+1)						      ! asdf
				   
	    n_c_tot = n_c_tot + 1									      ! Also, we have to add the Cell to our total (EXACT) Count.
	    
      end do												    ! (END "for every ROW of the SUBfile")

      10 close (606+j, status='delete')									    ! *: Since we couldn't know how many rows there were in the SUBfile beforehand,
      								!we couldn't specify a counter to determine the end of the "do" cycle; instead, we told the program "end=10", which means
      								! "go to '10' when you reach the end of the SUBfile. Well, '10' is THIS LINE, so when we reach the end of the SUBfile during
      								! the "do" cycle (which means that we're done reading it), we delete the SUBfile and go on
    
    end do												  ! (END "for every SUBsnapshot")
  
    ! Now that we've read (and deleted) every SUBsnapshot -and consequently built the whole "Single Vectors"-, we go ahead
    ! and write the FULL Snapshot, starting by ordering the Single Vectors by the Cells' distance to origin and "z" coordinate:
    
    allocate(Snap % puesto(n_c_tot)); Snap % puesto = 0.
    allocate(Snap % equis2(n_c_tot)); Snap % equis2 = 0.
    allocate(Snap % ygriega2(n_c_tot)); Snap % ygriega2 = 0.
    allocate(Snap % ceta2(n_c_tot)); Snap % ceta2 = 0.

    Snap % equis2(:) = HUGE
    Snap % ygriega2(:) = HUGE
    
    do i = 1, n_c_tot
      Snap % puesto(i) = i
      Snap % equis2(i) = min(Snap % equis2(i),  &
                     Snap % x_snap(i))
      Snap % ygriega2(i) = min(Snap % ygriega2(i),  &
                     Snap % y_snap(i))
      Snap % ceta2(i) = Snap % z_snap(i)
    end do
    
  call Sort % Three_Real_Carry_Int(Snap % equis2, Snap % ygriega2, Snap % ceta2, Snap % puesto)		! <------(This is the function that will do the sorting for us)										
    
    select case(ts)											  ! We then set the NAME of the SnapShot, with the same syntax as
    													  ! the vtu's (because it looks cool). To get that, we condition
    													  ! the format depending on the OoM of the TS at which we are:
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

    open(unit=607,file = trim(filename),form='formatted',status='unknown')				  ! So yet another time, let's go ahead and create the file...
    
    do i = 1, n_c_tot											  ! ...and, for every Cell in the WHOLE Domain...
      
      write (607, *) Snap % x_snap(Snap % puesto(i)), ',',     &					    ! ...write a row of the file containing its "x" Coordinate,...
        	     Snap % y_snap(Snap % puesto(i)), ',',     &					    ! ..."y" Coordinate,...
        	     Snap % z_snap(Snap % puesto(i)), ',',     &					    ! ..."z" Coordinate,...
          	     Snap % u_snap(Snap % puesto(i)), ',',     &					    ! ...Streamwise Velocity,...
        	     Snap % v_snap(Snap % puesto(i)), ',',     &					    ! ...Vertical Velocity,...
        	     Snap % w_snap(Snap % puesto(i)), ',',     &					    ! ...Spanwise Velocity,...
        	     Snap % p_snap(Snap % puesto(i)), ',',     &					    ! ...and Pressure:
		     Snap % dudx_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dudy_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dudz_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dvdx_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dvdy_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dvdz_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dwdx_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dwdy_snap(Snap % puesto(i)), ',',  &					    ! asdf
		     Snap % dwdz_snap(Snap % puesto(i))							    ! asdf
        								     
    end do												  ! (END "for every Cell in the WHOLE Domain")
    
    deallocate(Snap % u_snap)										  ! And we're pretty much done; all there's left is to deallocate the...
    deallocate(Snap % v_snap)										  ! ...space that we reserved for the "Single Vectors" (because we don't...
    deallocate(Snap % w_snap)										  ! ...need them anymore, since we've already saved the corresponding...
    deallocate(Snap % p_snap)										  ! ...information in the SNAPSHOT),...
    deallocate(Snap % x_snap)										  ! -----------------------------------
    deallocate(Snap % y_snap)										  ! -----------------------------------
    deallocate(Snap % z_snap)										  ! -----------------------------------
    deallocate(Snap % puesto)										  ! -----------------------------------
    deallocate(Snap % equis2)										  ! -----------------------------------
    deallocate(Snap % ygriega2)										  ! -----------------------------------
    deallocate(Snap % ceta2)										  !  ------------------------------------------------------------------
    deallocate(Snap % dudx_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dudy_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dudz_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dvdx_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dvdy_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dvdz_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dwdx_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dwdy_snap)									  !  ------------------------------------------------------------------
    deallocate(Snap % dwdz_snap)									  !  --------------------------------------------------------------------------
      
    close(607)												  ! ...close the SnapShot File,...

    n_c_tot = 0;											  ! ...and, FINALLY reset the Exact Cell Counter to 0 (not sure if it's needed)

  end if												! [END "if we HAVE started to (...) in a SINGLE SnapShot"]
  
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ENDING YVES ZONE (2/2)



    !--------------!
    !   Velocity   !
    !--------------!
    call Results % Save_Vector_Real("Velocity [m/s]", plot_inside,  &
                                    Flow % u % n(c_f:c_l),          &
                                    Flow % v % n(c_f:c_l),          &
                                    Flow % w % n(c_f:c_l),          &
                                    f8, f9, data_offset, run)

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
        call Results % Save_Scalar_Real("Turbulent Heat Flux X [K m/s]",     &
                                        plot_inside,                         &
                                        Turb % ut % n(c_f:c_l),              &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Turbulent Heat Flux Y [K m/s]",     &
                                        plot_inside,                         &
                                        Turb % vt % n(c_f:c_l),              &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Turbulent Heat Flux Z [K m/s]",     &
                                        plot_inside,                         &
                                        Turb % wt % n(c_f:c_l),              &
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
                                      plot_inside,                          &
                                      kin_vis_t(c_f:c_l),                   &
                                      f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Results % Save_Scalar_Real("Reynolds Stress XX [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % uu % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Reynolds Stress YY [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % vv % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Reynolds Stress ZZ [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % ww % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Reynolds Stress XY [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % uv % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Reynolds Stress XZ [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % uw % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Reynolds Stress YZ [m^2/s^2]",  &
                                      plot_inside,                     &
                                      Turb % vw % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("Turbulent Heat Flux X [K m/s]",  &
                                        plot_inside,                      &
                                        Turb % ut % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Turbulent Heat Flux Y [K m/s]",  &
                                        plot_inside,                      &
                                        Turb % vt % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Turbulent Heat Flux Z [K m/s]",  &
                                        plot_inside,                      &
                                        Turb % wt % n(c_f:c_l),           &
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
      do c1 = c_f, c_l
        save_01(c1) = Turb % uu_res(c1) - Turb % u_mean(c1) * Turb % u_mean(c1)
        save_02(c1) = Turb % vv_res(c1) - Turb % v_mean(c1) * Turb % v_mean(c1)
        save_03(c1) = Turb % ww_res(c1) - Turb % w_mean(c1) * Turb % w_mean(c1)
        save_04(c1) = Turb % uv_res(c1) - Turb % u_mean(c1) * Turb % v_mean(c1)
        save_05(c1) = Turb % uw_res(c1) - Turb % u_mean(c1) * Turb % w_mean(c1)
        save_06(c1) = Turb % vw_res(c1) - Turb % v_mean(c1) * Turb % w_mean(c1)
      end do
      call Results % Save_Scalar_Real("Mean Reynolds Stress XX [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_01(c_f:c_l),                     &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Mean Reynolds Stress YY [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_02(c_f:c_l),                     &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Mean Reynolds Stress ZZ [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_03(c_f:c_l),                     &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Mean Reynolds Stress XY [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_04(c_f:c_l),                     &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Mean Reynolds Stress XZ [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_05(c_f:c_l),                     &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("Mean Reynolds Stress YZ [m^s/s^2]",  &
                                      plot_inside,                          &
                                      save_06(c_f:c_l),                     &
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
        call Results % Save_Scalar_Real("Mean Turbulent Heat Flux X [K m/s]",  &
                                        plot_inside,                           &
                                        save_01(c_f:c_l),                      &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Mean Turbulent Heat Flux Y [K m/s]",  &
                                        plot_inside,                           &
                                        save_02(c_f:c_l),                      &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("Mean Turbulent Heat Flux Z [K m/s]",  &
                                        plot_inside,                           &
                                        save_03(c_f:c_l),                      &
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

    ! Save y+ for all turbulence models
    if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       Turb % model .ne. DNS) then
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

  call Cpu_Timer % Stop('Save_Vtu_Results')

  end subroutine
