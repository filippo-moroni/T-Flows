!==============================================================================!
  subroutine Control_Mod_Max_Gauss_Gradients_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_GAUSS_GRADIENTS_ITERATIONS',  &
                                  val, val, verbose)

  end subroutine