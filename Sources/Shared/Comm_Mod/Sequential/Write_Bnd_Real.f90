!==============================================================================!
  subroutine Write_Bnd_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" boundary-cell-based array    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer(DP)      :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Write "distributed" boundary cell data 
  do c = 1, Comm % nb_tot
    write(fh) array(c)
  end do

  disp = disp + Comm % nb_tot * RP

  end subroutine
