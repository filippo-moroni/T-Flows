!==============================================================================!
  subroutine Finalize_Work(Work)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
!==============================================================================!

  print '(a)', ' #===========================#'
  print '(a)', ' #       Work_Mod usage      #'
  print '(a)', ' #---------------------------#'
  print '(a,i2,a3,i2,a2)', ' # Real cell arrays: ', Work % max_r_cell, ' / ',  &
                                                    Work % req_r_cell, ' #'
  print '(a,i2,a3,i2,a2)', ' # Real face arrays: ', Work % max_r_face, ' / ',  &
                                                    Work % req_r_face, ' #'
  print '(a,i2,a3,i2,a2)', ' # Real node arrays: ', Work % max_r_node, ' / ',  &
                                                    Work % req_r_node, ' #'
  print '(a,i2,a3,i2,a2)', ' # Int. cell arrays: ', Work % max_i_cell, ' / ',  &
                                                    Work % req_i_cell, ' #'
  print '(a,i2,a3,i2,a2)', ' # Int. face arrays: ', Work % max_i_face, ' / ',  &
                                                    Work % req_i_face, ' #'
  print '(a,i2,a3,i2,a2)', ' # Int. node arrays: ', Work % max_i_node, ' / ',  &
                                                    Work % req_i_node, ' #'
  print '(a)', ' #---------------------------#'

  end subroutine

