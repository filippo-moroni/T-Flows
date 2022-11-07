!==============================================================================!
  subroutine Create_Icosahedron(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, is
  real    :: a, t
!==============================================================================!

  t = (1.0+(5.0)**0.5)/2.0
  a = (1.0+t**2.0)**0.5

  Pol % n_faces = 20
  Pol % n_nodes = 12

  i = 1
  Pol % nodes_xyz(i,1) = t/a
  Pol % nodes_xyz(i,2) = 1.0/a
  Pol % nodes_xyz(i,3) = 0.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = -t/a
  Pol % nodes_xyz(i,2) = 1.0/a
  Pol % nodes_xyz(i,3) = 0.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = t/a
  Pol % nodes_xyz(i,2) = -1.0/a
  Pol % nodes_xyz(i,3) = 0.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = -t/a
  Pol % nodes_xyz(i,2) = -1.0/a
  Pol % nodes_xyz(i,3) = 0.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = 1.0/a
  Pol % nodes_xyz(i,2) = 0.0/a
  Pol % nodes_xyz(i,3) = t/a
  i = i+1
  Pol % nodes_xyz(i,1) = 1.0/a
  Pol % nodes_xyz(i,2) = 0.0/a
  Pol % nodes_xyz(i,3) = -t/a
  i = i+1
  Pol % nodes_xyz(i,1) = -1.0/a
  Pol % nodes_xyz(i,2) = 0.0/a
  Pol % nodes_xyz(i,3) = t/a
  i = i+1
  Pol % nodes_xyz(i,1) = -1.0/a
  Pol % nodes_xyz(i,2) = 0.0/a
  Pol % nodes_xyz(i,3) = -t/a
  i = i+1
  Pol % nodes_xyz(i,1) = 0.0/a
  Pol % nodes_xyz(i,2) = t/a
  Pol % nodes_xyz(i,3) = 1.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = 0.0/a
  Pol % nodes_xyz(i,2) = -t/a
  Pol % nodes_xyz(i,3) = 1.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = 0.0/a
  Pol % nodes_xyz(i,2) = t/a
  Pol % nodes_xyz(i,3) = -1.0/a
  i = i+1
  Pol % nodes_xyz(i,1) = 0.0/a
  Pol % nodes_xyz(i,2) = -t/a
  Pol % nodes_xyz(i,3) = -1.0/a
  do is = 1,Pol % n_faces
     Pol % faces_n_nodes(is) = 3
  end do
  is = 1
  i = 1
  Pol % faces_n(is,i) = 1
  i = i+1
  Pol % faces_n(is,i) = 9
  i = i+1
  Pol % faces_n(is,i) = 5
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 1
  i = i+1
  Pol % faces_n(is,i) = 6
  i = i+1
  Pol % faces_n(is,i) = 11
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 3
  i = i+1
  Pol % faces_n(is,i) = 5
  i = i+1
  Pol % faces_n(is,i) = 10
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 3
  i = i+1
  Pol % faces_n(is,i) = 12
  i = i+1
  Pol % faces_n(is,i) = 6
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 2
  i = i+1
  Pol % faces_n(is,i) = 7
  i = i+1
  Pol % faces_n(is,i) = 9
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 2
  i = i+1
  Pol % faces_n(is,i) = 11
  i = i+1
  Pol % faces_n(is,i) = 8
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 4
  i = i+1
  Pol % faces_n(is,i) = 10
  i = i+1
  Pol % faces_n(is,i) = 7
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 4
  i = i+1
  Pol % faces_n(is,i) = 8
  i = i+1
  Pol % faces_n(is,i) = 12
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 1
  i = i+1
  Pol % faces_n(is,i) = 11
  i = i+1
  Pol % faces_n(is,i) = 9
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 2
  i = i+1
  Pol % faces_n(is,i) = 9
  i = i+1
  Pol % faces_n(is,i) = 11
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 3
  i = i+1
  Pol % faces_n(is,i) = 10
  i = i+1
  Pol % faces_n(is,i) = 12
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 4
  i = i+1
  Pol % faces_n(is,i) = 12
  i = i+1
  Pol % faces_n(is,i) = 10
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 5
  i = i+1
  Pol % faces_n(is,i) = 3
  i = i+1
  Pol % faces_n(is,i) = 1
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 6
  i = i+1
  Pol % faces_n(is,i) = 1
  i = i+1
  Pol % faces_n(is,i) = 3
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 7
  i = i+1
  Pol % faces_n(is,i) = 2
  i = i+1
  Pol % faces_n(is,i) = 4
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 8
  i = i+1
  Pol % faces_n(is,i) = 4
  i = i+1
  Pol % faces_n(is,i) = 2
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 9
  i = i+1
  Pol % faces_n(is,i) = 7
  i = i+1
  Pol % faces_n(is,i) = 5
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 10
  i = i+1
  Pol % faces_n(is,i) = 5
  i = i+1
  Pol % faces_n(is,i) = 7
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 11
  i = i+1
  Pol % faces_n(is,i) = 6
  i = i+1
  Pol % faces_n(is,i) = 8
  is = is+1
  i = 1
  Pol % faces_n(is,i) = 12
  i = i+1
  Pol % faces_n(is,i) = 8
  i = i+1
  Pol % faces_n(is,i) = 6

  end