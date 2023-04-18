!==============================================================================!
  subroutine Calculate_Inertia_Prime(Grid)
!------------------------------------------------------------------------------!
! Calculates the inertia prime tensor for the TVM model.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
!-----------------------------------[Locals]-----------------------------------!
    
    real    :: half_t    ! Half trace of the inertia tensor
    integer :: c
        
!==============================================================================!

  ! Perform the calculation through all the cells
  do c = 1, Grid % n_cells
  
  half_t = 0.5*(Grid % ixx(c) + Grid % iyy(c) + Grid % izz(c))   ! Half of the trace of the inertia tensor
  
  Grid % ixxp(c) = Grid % ixx(c) - half_t                        ! The inertia prime tensor has only modified diagonal components
  Grid % iyyp(c) = Grid % iyy(c) - half_t
  Grid % izzp(c) = Grid % izz(c) - half_t
  
  Grid % ixyp(c) = Grid % ixy(c)                                 ! Deviatoric components are the same
  Grid % ixzp(c) = Grid % ixz(c)
  Grid % iyzp(c) = Grid % iyz(c)
                                      
  end do

  end subroutine
