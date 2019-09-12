
!
! modules:     cnstr_data:
!                    sub:  set_space: set_overlap_file(1), set_config(2), set_cnstr_array(3)
!                ite_data:
!                    sub: set_ite_data   
!

MODULE cnstr_data

  USE PRECISION,only:dpk
  IMPLICIT NONE

  PUBLIC
  CHARACTER(len=16), ALLOCATABLE, DIMENSION(:) :: overlap
  INTEGER,           ALLOCATABLE, DIMENSION(:) :: ix, iy , ndi, ic, idcs
  INTEGER,           ALLOCATABLE, DIMENSION(:) :: nhf, lhf, ll, nmin, nmax, is, noll
  REAL(dpk),       ALLOCATABLE, DIMENSION(:,:) :: cvec

  INTERFACE set_space
     MODULE PROCEDURE set_overlap_file, set_cnstr_array, set_config
  END INTERFACE

  PRIVATE set_overlap_file, set_cnstr_array, set_config
CONTAINS 
  !S
  !S
  !S
  SUBROUTINE set_overlap_file(nl)
    !
    IMPLICIT NONE
    INTEGER  nl
    !
    ALLOCATE(overlap(nl))
    !
    WRITE(*,*) '# h2e::set_overlap_file:: allocation for overlap matrices done. nl = ', nl

  END SUBROUTINE set_overlap_file
  !S
  !S
  !S
  SUBROUTINE set_config(nd)
    !
    IMPLICIT NONE
    INTEGER  nchl, nd
    !
    ALLOCATE( ndi(nd)  )
    ALLOCATE( nhf(nd)  )
    ALLOCATE( lhf(nd)  )
    ALLOCATE( ll(nd)   )
    ALLOCATE( nmin(nd) )
    ALLOCATE( nmax(nd) )
    ALLOCATE( noll(nd) )
    ALLOCATE( is(nd)   ) 
    ALLOCATE( idcs(nd) )
    !
    WRITE(*,'(a40,i10)') '# h2e::set_config: set configurations done.'
    WRITE(*,'(a40,i10)') '# h2e::                               nd = ', nd 
    
  END SUBROUTINE set_config
  !S
  !S
  !S
  SUBROUTINE set_cnstr_array(nocnstr, nmx)
    !
    IMPLICIT NONE
    INTEGER nocnstr, nmx
    !
    ALLOCATE(    ix(nocnstr)      ) 
    ALLOCATE(    iy(nocnstr)      ) 
    ALLOCATE(    ic(nocnstr)      ) 
    ALLOCATE(  cvec(nmx, nocnstr) )
    !
    WRITE(*,'(a40,i10)') '# h2e::set_cnstr_array: set constraints array done.'
    WRITE(*,'(a40,i10)') '# h2e::                                   noc = ', nocnstr 
    !
  END SUBROUTINE set_cnstr_array

END MODULE cnstr_data
!M
!M
!M
MODULE ite_data
  !
  USE PRECISION, ONLY:dpk
  !
  IMPLICIT NONE
  PUBLIC
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: cx, bx
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: er, vec 
  REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: finit

  INTERFACE ite_space
     MODULE PROCEDURE set_ite_data
  END INTERFACE

  PRIVATE set_ite_data

CONTAINS


  !S
  !S
  !S



  SUBROUTINE set_ite_data(ndim, nchl)
    IMPLICIT NONE
    INTEGER  ndim, nchl
!!!....................

    WRITE(*,'(a40,i10)') '# h2e::set_ite:  allocation for energy matrix done. '
    WRITE(*,'(a40,i10)') '# h2e::set_ite:                              ndim = ', ndim 
    WRITE(*,'(a40,i10)') '# h2e::set_ite:                              nchl = ', nchl 

    ALLOCATE(er(ndim,nchl))
    ALLOCATE(vec(ndim,nchl))
    ALLOCATE(finit(ndim))

  END SUBROUTINE set_ite_data


END MODULE ite_data
!eof
