!==============================================================================!
  subroutine Save_Scalar_Real(grid, var_name, plot_inside,  &
                              val, fs, fp,                  &
                              data_offset, sweep)
!------------------------------------------------------------------------------!
!   Writes one real scalar defined over cells.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: var_name
  logical          :: plot_inside     ! plot results inside?
  real             :: val(-grid % n_bnd_cells:grid % n_cells)
  integer          :: fs, fp          ! file unit sequential and parallel
  integer          :: data_offset
  integer          :: sweep           ! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(4)         :: data_size
  integer            :: c, c2
  character(len=80)  :: str1
  integer, parameter :: IP=8, RP=8, SP=4
!==============================================================================!

  data_size = 0

  ! Header
  if(sweep .eq. 1) then
    if(n_proc > 1 .and. this_proc .eq. 1) then
      write(fs) IN_4                                  //  &
                '<PDataArray type="Float64"'          //  &
                ' Name="' // trim(var_name) // '"/>'  // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type="Float64"'        //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = grid % n_cells * RP
      write(fp) data_size
      do c = 1, grid % n_cells
        write(fp) val(c)
      end do
    else
      do c2 = -grid % n_bnd_cells, -1
        data_size = data_size + RP
      end do
      write(fp) data_size
      do c2 = -grid % n_bnd_cells, -1
        write(fp) val(c2)
      end do
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset  &
                  + grid % n_cells * RP
    else
      do c2 = -grid % n_bnd_cells, -1
        data_offset = data_offset + RP
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine