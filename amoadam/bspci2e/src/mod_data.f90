!MOD!
!
!
!
!
!
MODULE DATA

  USE PRECISION, ONLY:dpk

  PUBLIC write_en, write_v_mx, read_v_mx, write_wf1e, read_v, read_mx 


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! store matrix elements ( upper - diagonal )
!
  SUBROUTINE write_mx(l, mx, filename, mx_type)
    !
    !    USE PRECISION,  ONLY : dpk
    USE UTILS,      ONLY : datafile
    !
    IMPLICIT NONE
    !arg!
    INTEGER                   :: l
    REAL(DPK), DIMENSION(:,:) :: mx
    CHARACTER(LEN=*)          :: filename
    CHARACTER(len=6),optional :: mx_type
    !loc!
    INTEGER                   :: ifile
    INTEGER                   :: dim_1, dim_2 
    INTEGER                   :: I, J
    !exe!

    WRITE(*,'(a60)') " subroutine write_mx in."
    !

    ifile= 56
    
    !,,,,
    
    dim_1 = SIZE(mx,1)
    dim_2 = SIZE(mx,2)
    CALL datafile(ifile, l, filename)
     
    WRITE(ifile) dim_1, dim_2    
    
    IF(PRESENT(mx_type)) THEN
       DO I = 1, dim_1        
          WRITE(ifile) (mx(i,j), j = 1, dim_2 )
       ENDDO
       WRITE(*, '(a60,a10)') "mx_type = ", mx_type
    ELSE
       DO i = 1, dim_1        
          WRITE(ifile) (mx(i,j), j = i, dim_2 )
!          IF(i==10) THEN
!             WRITE(*,'(10f10.3)') (mx(i,j),j=i,dim_2)
!          ENDIF          
       ENDDO
    ENDIF        
    CLOSE(ifile)

    WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
    WRITE(*,'(a60)') "subroutine write_mx out.", filename
    RETURN
    
  END SUBROUTINE WRITE_MX
!
! read matrix elements
!
  SUBROUTINE read_mx(l, mx, filename, mx_type)
    !mod!
    !    USE PRECISION,  ONLY : dpk
    USE UTILS,      ONLY : datafile, check
    !
    IMPLICIT NONE
    !arg!
    INTEGER                   :: l
    REAL(DPK), DIMENSION(:,:) :: mx
    CHARACTER(LEN=*)          :: filename
    CHARACTER(len=6),OPTIONAL :: mx_type
    !loc!
    INTEGER                   :: ifile
    INTEGER                   :: dim_1, dim_2, ndim_1, ndim_2 
    INTEGER                   :: i, j
    INTEGER                   :: ierror
    !exe!
    WRITE(*,'(a60)') " subroutine read_mx in."
    WRITE(*,'(a60)') " filename = ", filename
    
    ierror = 0
    ifile  = 56

    
    !,,,,,,
    
    ndim_1 = SIZE(mx,1)
    ndim_2 = SIZE(mx,2)

    CALL DATAFILE(ifile, l, filename)
     
    READ(ifile) dim_1, dim_2
    
    call check(ndim_1, dim_1, ierror)
    call check(ndim_2, dim_2, ierror+1)

    IF(PRESENT(mx_type)) THEN
       store_full:DO i = 1, dim_1
          READ(ifile) (mx(i,j), j = 1, dim_2 )
       ENDDO store_full
       WRITE(*, '(a60,a10)') "mx_type = ", mx_type
    ELSE
       store_upper:DO i = 1, dim_1
          READ(ifile) (mx(i,j), j = i, dim_2 )
       ENDDO store_upper
    ENDIF

    CLOSE(ifile)

    WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
    WRITE(*,'(a60)') " subroutine read_mx out."
    RETURN
    
  END SUBROUTINE read_mx
  !
  ! pointer_version
  !
  SUBROUTINE READ_MX_P(l, mx, filename, mx_type)

    !    USE PRECISION,  ONLY : dpk
    USE UTILS,      ONLY : datafile
    !
    IMPLICIT NONE
    !.... Arguments
    INTEGER                            :: l
    REAL(DPK), DIMENSION(:,:), POINTER :: mx
    CHARACTER(LEN=*)                   :: filename
    CHARACTER(len=6),OPTIONAL          :: mx_type
    !loc!
    INTEGER                            :: ifile
    INTEGER                            :: i,j
    INTEGER                            :: dim_1, dim_2 
    !exe!

    WRITE(*,'(a60)') " subroutine read_mx_p in."
    WRITE(*,'(a60)') " filename = ", filename
    
    
    ifile = 56

    
    !........................
    
    CALL DATAFILE(ifile, l, filename)
     
    READ(ifile) dim_1, dim_2

    ALLOCATE(mx(dim_1,dim_2))

    IF(PRESENT(mx_type)) THEN
       read_full:DO I = 1, dim_1
          READ(ifile) (mx(i,j), j = 1, dim_1 )
       ENDDO read_full
       WRITE(*, '(a60,a10)') "mx_type = ", mx_type
    ELSE
       read_upper:DO I = 1, dim_1
          READ(ifile) (mx(i,j), j = i, dim_1 )
       ENDDO read_upper
    ENDIF    
    CLOSE(ifile)
    WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
    WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
    WRITE(*,'(a60)') " subroutine read_mx_p out."
    RETURN
    
  END SUBROUTINE READ_MX_P

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! store vectors
SUBROUTINE WRITE_V(L, V, FILENAME)

  !  use precision,  only : DPK
  use UTILS,      only : ASCIIFILE
  use PARAM,      only : nb, nbs, ncs


  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:) :: V
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I, J
  integer dim 
  integer er
  !.........................


  FILE = 56

!........................

  dim = size(V)
  
  CALL ASCIIFILE(FILE, L, FILENAME)
     
  WRITE(FILE,'(i5)') dim

  DO I = 1, dim
     WRITE(FILE,'(i5,2x,E25.14)') i, V(i) 
  ENDDO
     
  CLOSE(FILE)
  RETURN

END SUBROUTINE WRITE_V
!
!
!
SUBROUTINE READ_V(L, V, FILENAME)

  !  use precision,  only : DPK
  use UTILS,      only : ASCIIFILE
  use PARAM,      only : nb, nbs, ncs


  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:) :: V
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I, J
  integer dim 
  integer er
  !.........................


  FILE = 56

!........................

  dim = size(V)
  
  CALL ASCIIFILE(FILE, L, FILENAME)
     
  READ(FILE,'(i5)') dim

  DO i = 1, dim
     READ(FILE,'(i5,2x,E25.14)') j, V(i) 
  ENDDO
     
  CLOSE(FILE)
  RETURN
  
END SUBROUTINE READ_V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! write/read  E/C(E) for each l = 0,1,2
!
!
!SUBROUTINE write_v_mx(l, v, mx, filename, mx_type)

  SUBROUTINE write_v_mx(l, v, mx, filename)
  !mod!
    !  USE PRECISION,  ONLY: dpk
  USE UTILS,      ONLY: datafile, check
  !
  implicit none
  !arg!
  integer L
  real(DPK), dimension(:)   :: v
  real(DPK), dimension(:,:) :: mx
  character(LEN=*)          :: filename
!  CHARACTER(len=6),OPTIONAL :: mx_type
  !loc!
  integer                   :: dim_v, dim_1, dim_2 
  !
  integer                   :: ifile
  integer                   :: ierror
  integer                   :: i, j
  !exe!
  WRITE(*,'(a60)') " subroutine write_v_mx in."
  ierror =  1
  ifile  = 56

  !,,, set dims

  dim_v = SIZE(v)
  dim_1 = SIZE(mx,1)
  dim_2 = SIZE(mx,2)

  !
  CALL check(dim_v, dim_1, ierror)
  CALL check(dim_v, dim_2, ierror+1)
  !
 
  CALL datafile(ifile, l, filename)     
  WRITE(ifile) dim_1, dim_2  
  !  IF(PRESENT(mx_type)) THEN
  store_v_mx_full: DO i = 1, dim_1
     WRITE(ifile)  v(i)  
     WRITE(ifile)( mx(i,j), j = 1, dim_2 )
  ENDDO store_v_mx_full
  
  !  ELSE
  !     store_v_mx_upper: DO i = 1, dim_1
  !        WRITE(ifile)  v(i)  
  !        WRITE(ifile)( mx(i,j), j = i, dim_2 )
  !     ENDDO store_v_mx_upper
  
  !  ENDIF
  

  CLOSE(ifile)

  WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
  WRITE(*,'(a60)') " subroutine write_v_mx out.", filename
  RETURN
  
END SUBROUTINE write_v_mx
!.......................
!SUBROUTINE read_v_mx(l, v, mx, filename, my_type)
  SUBROUTINE read_v_mx(l, v, mx, filename)
  !mod!
    !  USE PRECISION,  ONLY : dpk
  USE UTILS,      ONLY : datafile
  !
  IMPLICIT NONE
  !arg!
  INTEGER L
  REAL(DPK), DIMENSION(:),   POINTER :: v
  REAL(DPK), DIMENSION(:,:), POINTER :: mx
  CHARACTER(LEN=*)                   :: filename
!  CHARACTER(len=6),OPTIONAL          :: mx_type
  !loc!
  INTEGER                            :: ifile
  INTEGER                            :: dim_1, dim_2 
  INTEGER                            :: ierror
  INTEGER                            :: i,j
  !exe!
  !,,,,,,,
  WRITE(*,'(a60)') " subroutine read_v_mx in."
  WRITE(*,'(a60)') " filename = ", filename
  ierror = 1
  ifile = 56

!,,,,,

 
  CALL datafile(ifile, l, filename)

  READ(ifile) dim_1, dim_2
  !     
  ALLOCATE( v( dim_1 ) )
  ALLOCATE( mx( dim_1, dim_2 ) )
  !  
  v  = 0.0_dpk
  mx = 0.0_dpk

  !  IF(PRESENT(mx_type)) THEN
  read_v_mx_full: DO i = 1, dim_1
     READ(ifile)  v(i)  
     READ(ifile)( mx(i,j), j = 1, dim_2 )
  ENDDO read_v_mx_full
     !  ELSE
     !      fac = 1     
     !     read_v_mx_upper: DO i = 1, dim_1
     !        read(ifile)  v(i)  
     !        read(ifile)( mx(i,j), j = i, dim_2 )
     !   mx(j,i) = fac * mx(i,j)
     !     ENDDO read_v_mx_upper
     
     !  ENDIF


  CLOSE(ifile)

  WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
  WRITE(*,'(a60)') " subroutine read_v_mx out."
  RETURN
END SUBROUTINE read_v_mx


END MODULE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!EOF
