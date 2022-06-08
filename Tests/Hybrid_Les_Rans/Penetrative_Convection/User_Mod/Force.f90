!==============================================================================!
  subroutine User_Mod_Force(Flow, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for velocity.      !
!   It is called from "Compute_Velocity" function, just before calling the     !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)   :: Flow
  type(Var_Type)     :: ui        ! velocity component
  type(Matrix_Type)  :: a_matrix  ! system matrix
  real, dimension(:) :: b_vector  ! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  !----------------------------------------------------!
  !                                                    !
  !   Set source depending on the velocity component   !
  !                                                    !
  !----------------------------------------------------!

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'U' ) then
    do c = 1, Grid % n_cells
      b_vector(c) = b_vector(c) - static_cell_fx(c)
    end do
  end if

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'V' ) then
    do c = 1, Grid % n_cells
      b_vector(c) = b_vector(c) - static_cell_fy(c)
    end do
  end if

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'W' ) then
    do c = 1, Grid % n_cells
      b_vector(c) = b_vector(c) - static_cell_fz(c) * Grid % vol(c)
    end do
  end if

  end subroutine
