!==============================================================================!
  function N_Facets(Stl)
!------------------------------------------------------------------------------!
!   Returns number of facets                                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type) :: Stl
  integer         :: N_Facets
!==============================================================================!

  N_Facets = Stl % n_bnd_cells

  end function
