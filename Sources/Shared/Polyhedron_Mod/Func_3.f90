!==============================================================================!
  function Func_3(Pol, x, y, z)
!------------------------------------------------------------------------------!
!   Orthocircle surface [http://paulbourke.net/geometry/orthocircle/]          !
!   centered at (1.25,1.25,1.25)                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  real                   :: Func_3
  real                   :: x, y, z
!==============================================================================!

  x=x-1.25
  y=y-1.25
  z=z-1.25
  Func_3=((x**2+y**2-1)**2+z**2)*((y**2+z**2-1)**2+x**2)*((z**2+  &
         x**2-1)**2+y**2)-(0.075)**2*(1+3*(x**2+y**2+z**2))

  end
