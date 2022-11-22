!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Front_Type), pointer :: Front
  integer                   :: c, e, g, l, s, c1, c2, i_ele
  real                      :: cond_1, cond_2
  type(Var_Type)            :: t_0, t_1  ! introduced to hold gradients in
                                         ! each of the phases, needed for
                                         ! estimation of mass transfer in cells
!==============================================================================!

  ! Take aliases
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  ! Distinguish between liquid and vapor (Yohei doesn't like)
  call Vof % Get_Gas_And_Liquid_Phase(g, l)

  !------------------------------------------------!
  !   Compute gradients of temperature, imposing   !
  !    saturation temperature at the interface     !
  !------------------------------------------------!
  call Vof % Calculate_Grad_Matrix_With_Front()
  call Vof % Grad_Variable_With_Front(Flow % t, Vof % t_sat)

  ! First extrapolate temperature gradients from phase 0 to 1
  call Vof % Extrapolate_Normal_To_Front(Flow, t_1, towards=1)

  ! Then extrapolate temperature gradients from phase 1 to 0
  call Vof % Extrapolate_Normal_To_Front(Flow, t_0, towards=0)

  !--------------------------------------------------!
  !   Compute heat flux at the interface cell-wise   !
  !--------------------------------------------------!
  do c = 1, Grid % n_cells
    e = Vof % Front % elem_in_cell(c)

    if(e .gt. 0) then
      Vof % m_dot(c) = (  t_0 % x(c) * Front % Elem(e) % sx           &
                        + t_0 % y(c) * Front % Elem(e) % sy           &
                        + t_0 % z(c) * Front % Elem(e) % sz)          &
                        * Vof % phase_cond(2)                         &
                     - (  t_1 % x(c) * Front % Elem(e) % sx           &
                        + t_1 % y(c) * Front % Elem(e) % sy           &
                        + t_1 % z(c) * Front % Elem(e) % sz)          &
                        * Vof % phase_cond(1)
      Vof % m_dot(c) = Vof % m_dot(c) / 2.26e+6
    end if
  end do

  end subroutine
