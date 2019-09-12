!
!
!MODULE atom_2e      : 2-electron atomic system
!MODULE my_types     : 1-electron-atom/molecule
!
!
!
!
!
!
!
!
MODULE my_types
  !
  USE PRECISION, ONLY: dpk
  !
  IMPLICIT NONE
  !

  ! 
  ! atomic state quantum numbers (aqn)
  !

  !
  !
  TYPE atomic_system
     CHARACTER(len=20)  :: name
     REAL(dpk)          :: znuc
     CHARACTER(len=20)  :: potential
     CHARACTER(len=6)   :: gauge
  END TYPE atomic_system
  !
  !  
  TYPE state_index
     INTEGER :: n         ! principal qn
     INTEGER :: l         ! angular   qn 
     INTEGER :: m         ! magnetic  qn
     INTEGER :: p         ! parity     
     INTEGER :: n_l        ! number of coupled angular 
  END TYPE state_index



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TYPE molecular_system
     INTEGER                              :: n_atoms
     CHARACTER(len=20)                    :: name
     REAL(dpk), ALLOCATABLE, dimension(:) :: znuc
     CHARACTER(len=20)                    :: potential
     CHARACTER(len=6)                     :: gauge
  END TYPE molecular_system
  !
  !
  !


  !
  !
  TYPE basis
     INTEGER                              :: n      ! basis size
     INTEGER                              :: k      ! basis order
     REAL(dpk), ALLOCATABLE, DIMENSION(:) :: t      ! knot sequence
  END TYPE basis



  !  TYPE state
  !     !
  !     REAL(dpk)                          ::  ek     ! energy
  !   REAL(dpk), ALLOCATABLE, DIMENSION(:) ::  cb     ! coefficients  
  !   REAL(dpk)                            ::  dk     ! phase shift
  !   REAL(dpk)                            ::  ak     ! box-normalization
  !   REAL(dpk)                            ::  wk     ! surface amplitude
  !
  !END TYPE state


CONTAINS

  !
  ! calls input from bs1e_parameter.f90
  !

  SUBROUTINE init_basis(this)
    !
    USE PARAM,        ONLY: nb, kb, idbsp, rmax   
    USE PARAM,        ONLY: input                 ! read
    USE one_e_matrix, ONLY: mkgrid 
    !
    IMPLICIT NONE
    !ARG!
    TYPE(basis), INTENT(inout):: this
    !EXE!

    CALL input

    this%n = nb
    this%k = kb
    ALLOCATE( this%t(1:nb+kb) )
    CALL mkgrid(this%t)
    !
  END SUBROUTINE init_basis
  !S
  !S
  !S
  SUBROUTINE init_atomic_system(this)
    !
    IMPLICIT NONE
    !ARG!
    TYPE(atomic_system), INTENT(inout):: this
    !EXE!
    this%name      ='h'
    this%znuc      = 1.0_dpk
    this%potential ='coulombic'
    this%gauge     = 'v'
    !
  END SUBROUTINE init_atomic_system


  SUBROUTINE init_molecular_system(this)
    !
    IMPLICIT NONE
    !ARG!
    TYPE(molecular_system), INTENT(inout):: this
    !EXE!



    this%n_atoms   = 2
    this%name      ='h2p'
    !
    ALLOCATE(this%znuc(this%n_atoms))
    this%znuc      = 1.0_dpk
    !
    this%potential ='coulombic'
    this%gauge     = 'v'
    !
  END SUBROUTINE init_molecular_system



 
END MODULE my_types
