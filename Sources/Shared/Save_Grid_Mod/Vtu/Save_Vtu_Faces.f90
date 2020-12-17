!==============================================================================!
  subroutine Save_Vtu_Faces(grid, plot_shadows)
!------------------------------------------------------------------------------!
!   Writes boundary condition .faces.vtu or shadow .shadow.vtu file.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  logical, optional :: plot_shadows  ! plot shadow faces
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c2, n, s, s_f, s_l, cell_offset, data_offset, n_conns, fu
  character(SL) :: name_out, ext, str1, str2
  real          :: mag
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP = DP  ! int. precision is double precision
  integer, parameter :: RP = DP  ! real precision is double precision
!==============================================================================!

  ! Starting and ending counters; file extension
  s_f = 1
  s_l = grid % n_faces
  ext = '.faces.vtu'

  ! Fix counters and file extension if you are plotting shadows
  if(present(plot_shadows)) then
    if(plot_shadows) then
      s_f = minval(grid % faces_s(1:grid % n_faces),  &
              mask=grid % faces_s(1:grid % n_faces) .ne. 0)
      s_l = maxval(grid % faces_s(1:grid % n_faces),  &
              mask=grid % faces_s(1:grid % n_faces) .ne. 0)
      ext = '.shadows.vtu'
    end if
  end if

  if(s_l - s_f < 1) then
    print *, '# NOTE: No shadow faces in this domain, nothing to plot!'
    return
  end if

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do s = s_f, s_l
    n_conns = n_conns + grid % faces_n_nodes(s)
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File_Mod_Set_Name(name_out, processor=this_proc, extension=trim(ext))
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>'             // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid"'  //  &
                    ' version="0.1"'                    //  &
                    ' byte_order="LittleEndian">'       // LF
  write(fu) IN_1 // '<UnstructuredGrid>' // LF
  write(str1, '(i0.0)') grid % n_nodes
  write(str2, '(i0.0)') (s_l-s_f+1)
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"'  //  &
                    ' NumberOfCells="'        // trim(str2) // '">' // LF
  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + grid % n_nodes * RP * 3  ! prepare for next

  !-----------!
  !   Faces   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! Faces' nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP             ! prepare for next

  ! Faces' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Faces' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  !----------------------!
  !   The end of faces   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !---------------!
  !   Face data   !
  !---------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  ! Boundary conditions
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="BoundaryConditions"'     //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Number of nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="GridNumberOfNodes"'      //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Surface vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridSurfaceVectors"'     //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Surface normals
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridSurfaceNormals"'     //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Connection vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridConnectionVectors"'  //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu) IN_3 // '</CellData>'         // LF
  write(fu) IN_2 // '</Piece>'            // LF
  write(fu) IN_1 // '</UnstructuredGrid>' // LF

  !-------------------!
  !                   !
  !   Appended data   !
  !                   !
  !-------------------!
  write(fu) IN_0 // '<AppendedData encoding="raw">' // LF
  write(fu) '_'

  !-----------!
  !   Nodes   !
  !-----------!
  data_size = grid % n_nodes * RP * 3
  write(fu) data_size
  do n = 1, grid % n_nodes
    write(fu) grid % xn(n), grid % yn(n), grid % zn(n)
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! Faces' nodes
  data_size = n_conns * IP
  write(fu) data_size
  do s = s_f, s_l
    n = grid % faces_n_nodes(s)
    write(fu) grid % faces_n(1:n,s)-1
  end do

  ! Faces' offsets
  data_size = (s_l-s_f+1) * IP
  write(fu) data_size
  cell_offset = 0
  do s = s_f, s_l
    cell_offset = cell_offset + grid % faces_n_nodes(s)
    write(fu) cell_offset
  end do

  ! Faces' types
  data_size = (s_l-s_f+1) * IP
  write(fu) data_size
  do s = s_f, s_l
    if(grid % faces_n_nodes(s) .eq. 4) then
      write(fu) VTK_QUAD
    else if(grid % faces_n_nodes(s) .eq. 3) then
      write(fu) VTK_TRIANGLE
    else
      write(fu) VTK_POLYGON
    end if
  end do

  ! Boundary conditions
  ! (Check c1 and c2 for shadow faces, seems to be something messed up)
  data_size = (s_l-s_f+1) * IP
  write(fu) data_size
  do s = s_f, s_l
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      write(fu) grid % bnd_cond % color(c2)
    else
      write(fu) 0
    end if
  end do

  ! Number of nodes
  data_size = (s_l-s_f+1) * IP
  write(fu) data_size
  do s = s_f, s_l
    write(fu) grid % faces_n_nodes(s)
  end do

  ! Surface vectors
  data_size = (s_l-s_f+1) * RP * 3
  write(fu) data_size
  do s = s_f, s_l
    write(fu) grid % sx(s), grid % sy(s), grid % sz(s)
  end do

  ! Surface normals
  data_size = (s_l-s_f+1) * RP * 3
  write(fu) data_size
  do s = s_f, s_l
    mag = sqrt(grid % sx(s)**2 + grid % sy(s)**2 + grid % sz(s)**2)
    write(fu) grid % sx(s) / mag, grid % sy(s) / mag, grid % sz(s) / mag
  end do

  ! Connection vectors
  data_size = (s_l-s_f+1) * RP * 3
  write(fu) data_size
  do s = s_f, s_l
    write(fu) grid % dx(s), grid % dy(s), grid % dz(s)
  end do

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  end subroutine
