#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

!==============================================================================!
  module Petsc_Mod
!------------------------------------------------------------------------------!
!   Module used for PETSc linear solvers.                                      !
!                                                                              !
!   Note: This module has all member procedures bundled in one file.  When I   !
!         included them from separate files, I experienced difficulties with   !
!         include files above and definition of PETSc types.  Go figure?       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Iso_C_Binding
  use PetscVec,   only: tVec
  use PetscMat,   only: tMat
  use PetscKSP,   only: tKSP, tPC
  use Native_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  logical, parameter :: PETSC_ACTIVE = .true.
  integer, parameter :: OUT_OF_ITS   = -3      ! KSP_DIVERGED_ITS

  !----------------!
  !   Petsc type   !
  !----------------!
  type Petsc_Type

    type(Grid_Type), pointer :: pnt_grid

    ! Petsc-related variables
    type(tMat)                  :: A         ! sparse matrix
    type(tVec)                  :: x         ! solution vector
    type(tVec)                  :: b         ! right hand side
    type(tKSP)                  :: ksp       ! linear solver context
    type(tPc)                   :: pc        ! preconditioner
    type(PetscInt)              :: m_lower   ! unknowns in this proc.
    type(PetscInt)              :: m_upper   ! total number of unknowns
    type(PetscInt)              :: miter     ! maximum number of iterations
    type(PetscInt)              :: niter     ! performed number of iterations
    type(PetscInt), allocatable :: d_nnz(:)  ! diagonal stencil width per cell
    type(PetscInt), allocatable :: o_nnz(:)  ! off-diag. stencil width per cell
    type(PetscErrorCode)        :: err
    type(PetscEnum)             :: reason

    ! Global cell numbering for PETSc, which is
    ! different from T-Flows' and stars from zero   <---= IMPORTANT
    integer, allocatable :: glo(:)

    contains
      procedure :: Create_Petsc
      procedure :: Solve_Petsc

  end type

  contains

!==============================================================================!
  subroutine Create_Petsc(Pet, Nat, Grid)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
!for debug:  use Work_Mod, only: pet_glo => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Native_Type)        :: Nat
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, k, start
  integer, allocatable :: all_lower_ms(:)
!==============================================================================!

  Pet % pnt_grid => Grid

  if(this_proc < 2) print *, '# Initializing PETSc.'

  ! Total number of unknowns and unknowns in this processor only
  Pet % m_upper = Grid % comm % nc_tot
  Pet % m_lower = Grid % n_cells - Grid % comm % n_buff_cells

  !----------------------------------------+
  !    Create global numbering for PETSc   !
  !----------------------------------------+------------------!
  !    This one has little to do with global numbering from   !
  !    T-Flows and is unique for each number of processors    !
  !-----------------------------------------------------------!

  ! Dimensions must spread from all boundary cells through all ...
  ! ... buffers cells to successfully use Grid % Exchange_Cells_Int
  allocate(Pet % glo(-Grid % n_bnd_cells:Grid % n_cells))
  Pet % glo(:) = 0

  if(n_proc < 2) then
    Pet % glo(1:Grid % n_cells) = Grid % Comm % cell_glo(1:Grid % n_cells) - 1
  else
    start = 1  ! first row
    allocate(all_lower_ms(n_proc));  ! allocate array for all m_lowers
    all_lower_ms(:) = 0              ! important to initialize to zero

    ! Distribute m_lowers among all processors
    all_lower_ms(this_proc) = Pet % m_lower
    call Comm_Mod_Global_Sum_Int_Array(n_proc, all_lower_ms)

    start = sum(all_lower_ms(1:this_proc)) - Pet % m_lower

    ! Distribute global numbers over other processors
    do i = 1, Pet % m_lower
      Pet % glo(i) = i + start - 1
    end do
    call Grid % Exchange_Cells_Int(Pet % glo)

    !for debug:    ! Check what you got
    !for debug:    pet_glo(:) = Pet % glo(:)
    !for debug:    call Grid % Save_Debug_Vtu("petsc-numbering",      &
    !for debug:                               scalar_cell = pet_glo,  &
    !for debug:                               scalar_name = "Petsc Numbers [1]")
  end if

  !----------------------!
  !   Initialize PETSc   !
  !----------------------!
  call C_Petsc_Initialize()

  !--------------------------!
  !    Create PETSc matrix   !
  !--------------------------!
  call C_Petsc_Mat_Create(Pet % A)

  !----------------------------!
  !    Set PETSc matrix size   !
  !----------------------------!
  call C_Petsc_Mat_Set_Sizes(Pet % A, Pet % m_lower, Pet % m_upper)

  !---------------------------------------------------------!
  !   Set PETSc matrix type to MATAIJ (and pray it works)   !
  !---------------------------------------------------------!
  call C_Petsc_Mat_Set_Type_To_Mat_Aij(Pet % A)

  !-----------------------------------!
  !   Pre-allocate the PETSc matrix   !
  !-----------------------------------!

  ! Allocate memory for array with number of non-zero entries per row
  ! for entries in this processor (d_nnz), and other processors (o_nnz)
  allocate(Pet % d_nnz(Pet % m_lower))
  allocate(Pet % o_nnz(Pet % m_lower))
  Pet % d_nnz(:) = 0
  Pet % o_nnz(:) = 0

  ! Find number of nonzeros (nnz) in this processor and in other processors
  do i = 1, Pet % m_lower
    do j = Nat % A % row(i), Nat % A % row(i+1)-1
      k = Nat % A % col(j)
      if(Grid % Comm % cell_proc(k) .eq. this_proc) then
        Pet % d_nnz(i) = Pet % d_nnz(i) + 1
      else
        Pet % o_nnz(i) = Pet % o_nnz(i) + 1
      end if
    end do
  end do

  ! This will call both MPI and Seq versions of preallocation
  call C_Petsc_Mat_Aij_Set_Preallocation(Pet % A,        &
                                         Pet % d_nnz,    &
                                         Pet % o_nnz)

  !--------------------------!
  !   Create PETSc vectors   !
  !--------------------------!
  call C_Petsc_Vec_Create_Mpi(Pet % x, Pet % m_lower, Pet % m_upper)
  call C_Petsc_Vec_Create_Mpi(Pet % b, Pet % m_lower, Pet % m_upper)

  !-------------------------!
  !   Create PETSc solver   !
  !-------------------------!
  call C_Petsc_Ksp_Create(Pet % ksp)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine

!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         norm)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  character(*)         :: solver          ! solver
  character(*)         :: prec            ! preconditioner
  character(SL)        :: prec_opts(MSI)  ! preconditioner options
  type(Matrix_Type)    :: A
  real                 :: x(-Pet % pnt_grid % n_bnd_cells :  &
                             Pet % pnt_grid % n_cells)
  real                 :: b( Pet % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol      ! tolerance
  real,    intent(out) :: fin_res  ! final residual
  real,    optional    :: norm     ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(PetscInt) :: i, j, k, l
  character(SL)  :: solvers        ! fortran string to store solver
  character(SL)  :: precs          ! fortran string to store preconditioner
!==============================================================================!

  !-----------------------------------------------------------!
  !   Fill up PETSc matrix with values from original matrix   !
  !-----------------------------------------------------------!
  do i = 1, Pet % m_lower
    do j = A % row(i), A % row(i+1)-1
      k = A % col(j)
      call C_Petsc_Mat_Set_Value(Pet % A,           &  ! matrix
                                 Pet % glo(i),      &  ! row
                                 Pet % glo(k),      &  ! column
                                 A % val(j))           ! matrix entry
    end do
  end do

  ! The following two calls are needed after the calls to MatSetValue
  call C_Petsc_Assemble_Mat(Pet % A)

  !---------------------!
  !   Fill up vectors   ! (Pet % glo starts from zero)
  !---------------------!
  do i = 1, Pet % m_lower
    call C_Petsc_Vec_Set_Value(Pet % x, Pet % glo(i), x(i))
    call C_Petsc_Vec_Set_Value(Pet % b, Pet % glo(i), b(i))
  end do

  ! The following two calls are needed after the calls to VecSetValue
  call C_Petsc_Assemble_Vec(Pet % x)
  call C_Petsc_Assemble_Vec(Pet % b)

  !-----------------------------------!
  !   Set solver and preconditioner   !
  !-----------------------------------!
  solvers = solver;   l = len_trim(solvers);  solvers(l+1:l+1) = c_null_char
  precs   = prec;     l = len_trim(precs);    precs  (l+1:l+1) = c_null_char
  call C_Petsc_Set_Solver_And_Preconditioner(Pet % ksp,  &  ! solver
                                             Pet % pc,   &  ! preconditioner
                                             Pet % A,    &
                                             solvers,    &
                                             precs)

  !------------------------------------!
  !   Process preconditioner options   !
  !------------------------------------!
  if(prec_opts(1) .ne. '') then
    i = 1
    do while(i < MSI .and. prec_opts(i)(1:1) .ne. '')

      ! Option is just a single word (followed by another option or end)
      if( prec_opts(i)(1:1) .eq. '-' .and. prec_opts(i+1)(1:1) .eq. '-' .or. &
          prec_opts(i)(1:1) .eq. '-' .and. prec_opts(i+1)(1:1) .eq. '') then
        !debug: print *, 'A:', trim(prec_opts(i))
        call C_Petsc_Options_Value(trim(prec_opts(i)), "")
        i = i + 1

      ! Option is followed by a switch
      else
        !debug: print *, 'B:', trim(prec_opts(i)), ' ', trim(prec_opts(i+1))
        call C_Petsc_Options_Value(trim(prec_opts(i)), trim(prec_opts(i+1)))
        i = i + 2

      end if
    end do
  end if

  !---------------------------!
  !   Set solver tolerances   !
  !---------------------------!
  Pet % miter = miter
  call C_Petsc_Ksp_Set_Tolerances(Pet % ksp,     &
                                  tol,           &  ! PetscReal rtol
                                  tol,           &  ! PetscReal abstol
                                  Pet % miter)
  !-----------!
  !   Solve   !
  !-----------!
  call C_Petsc_Ksp_Solve(Pet % ksp, Pet % b, Pet % x)

  ! Check if converged
  call C_Petsc_Ksp_Converged_Reason(Pet % ksp, Pet % reason)

  ! Fetch the performed number of iterations
  call C_Petsc_Ksp_Get_Iteration_Number(Pet % ksp, Pet % niter)
  niter = Pet % niter

  ! Fetch the performed number of iterations
  call C_Petsc_Ksp_Get_Residual_Norm(Pet % ksp, fin_res)

  !-----------------------------------------------!
  !   Copy the solution back to T-Flows' vector   ! (Pet % glo starts from zero)
  !-----------------------------------------------!
  if(Pet % reason > 0 .or.  &            ! converged
     Pet % reason .eq. OUT_OF_ITS) then  ! simply ran out of iterations
    do i = 1, Pet % m_lower
      call C_Petsc_Vec_Get_Values(Pet % x,       &
                                  1,             &
                                  Pet % glo(i),  &
                                  x(i))
    end do
  else
    ! if(this_proc < 2) print *, ' # Warning: linear system failed to converge!'
  end if

  end subroutine

  end module