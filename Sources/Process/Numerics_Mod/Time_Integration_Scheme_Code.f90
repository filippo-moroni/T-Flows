!==============================================================================!
  integer function Numerics_Mod_Time_Integration_Scheme_Code(scheme_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: scheme_name
!==============================================================================!

  select case(scheme_name)

    case('LINEAR')
      Numerics_Mod_Time_Integration_Scheme_Code = LINEAR
    case('PARABOLIC')
      Numerics_Mod_Time_Integration_Scheme_Code = PARABOLIC
    case('RUNGE_KUTTA_3')
      Numerics_Mod_Time_Integration_Scheme_Code = RUNGE_KUTTA_3

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown time-integration scheme: ',  &
                 trim(scheme_name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop

  end select

  end function
