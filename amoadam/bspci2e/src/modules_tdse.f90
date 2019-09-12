!
MODULE deriv
  USE PRECISION
  IMPLICIT NONE
  REAL(DPK), POINTER, DIMENSION(:) :: yderiv
END MODULE deriv
!
!...
!
MODULE PARAMETER
  !
  USE PRECISION
  !
  IMPLICIT NONE
  !
  PUBLIC
  !common to all tdse*.f90
  INTEGER                                  :: nmax
  INTEGER                                  :: lmax
  INTEGER                                  :: ltot
  INTEGER                                  :: ntot
  INTEGER                                  :: neq  
  INTEGER,   ALLOCATABLE, DIMENSION(:)     :: n
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: en
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: dzr
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: dag
  INTEGER,   ALLOCATABLE, DIMENSION(:)     :: nsum       
  !
  !pulse!
  REAL(dpk)                                :: e0
  REAL(dpk)                                :: omeg
  REAL(dpk)                                :: tau
  REAL(dpk)                                :: cepd
  CHARACTER(len=10)                        :: spulse
  CHARACTER(len=6)                         :: sgauge
  !fcn 
  INTEGER                                  :: neval
  !
  ! free boundary condition (fr)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: dzi
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)   :: ch_index_en
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)   :: noch_en

  !

CONTAINS
!
! free
!
  SUBROUTINE make_space(nmax, ltot) 
    
    INTEGER nmax, ltot
    !          neq = nmax * ltot * 2
    
    ALLOCATE( n(ltot), nsum(ltot)    )
    ALLOCATE( dzr(nmax,nmax,ltot-1)  )
    ALLOCATE( dzi(nmax,nmax,ltot-1)  )
    ALLOCATE( en(nmax, ltot)         )
    ALLOCATE( ch_index_en(nmax, ltot))
    ALLOCATE( noch_en(nmax, ltot)    )
    ALLOCATE( dag(nmax*ltot*2)       )
    
  END SUBROUTINE make_space
  !
  !
  !

END MODULE PARAMETER
!eof
