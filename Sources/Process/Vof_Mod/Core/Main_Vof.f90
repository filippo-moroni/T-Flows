!==============================================================================!
  subroutine Main_Vof(Vof, Flow, turb, Sol, curr_dt)
!------------------------------------------------------------------------------!
!   Initialize inteface tracking simulations                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),     target :: Vof
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Solver_Type)           :: Sol
  integer, intent(in)         :: curr_dt     ! time step
!==============================================================================!

  if(Vof % model .eq. VOLUME_OF_FLUID) then

    !----------------------------------!
    !   Created a front or a surface   !
    !----------------------------------!
    if(Vof % track_front) then
      call Vof % Front % Place_At_Var_Value(Vof % fun,  &
                                            Sol,        &
                                            0.5,        &
                                            .true.)  ! don't print messages
      call Vof % Front % Print_Front_Statistics()
!     call Vof % Front % Save_Front(curr_dt)
    end if

    !----------------------!
    !   Create a surface   !
    !----------------------!
    if(Vof % track_surface) then
      call Surf_Mod_Place_At_Var_Value(Vof % Surf,  &
                                       Vof % fun,   &
                                       Sol,         &
                                       0.5,         &
                                       .true.)  ! don't print messages
      call Vof % Surf % Calculate_Curvatures_From_Elems()
      call Vof % Surf % Print_Surf_Statistics()
    end if

    !--------------------------------!
    !   Advance vof function (fun)   !
    !--------------------------------!
    call Update_Boundary_Values(Flow, turb, Vof, 'MULTIPHASE')
    call Vof % Compute_Vof(Sol, Flow % dt, curr_dt)

    !------------------------------------------------!
    !   Prepare smooth variant of the vof function   !
    !    for computation of normals and curvature    !
    !------------------------------------------------!
    if(Vof % surface_tension > TINY) then
      call Vof % Smooth_For_Curvature_Csf()
      call Vof % Curvature_Csf()
    end if

    !---------------------------------------------------------------!
    !   Update properties for other conservation equations to use   !
    !    (Maybe redundant since it is called from Main_Pro too.)    !
    !---------------------------------------------------------------!
    call Vof % Update_Physical_Properties()

  end if

  end subroutine
