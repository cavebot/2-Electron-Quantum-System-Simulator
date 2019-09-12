!!!
MODULE utils
  USE PRECISION, ONLY:DPK
  IMPLICIT NONE


  PUBLIC asciifile, datafile,check, print_mx, banded_mat_v_mul, dfile, afile, tafile
  
contains
! print out matrix
  SUBROUTINE print_mx(nmx, mx, smx, sf) 
    !
    USE PRECISION, ONLY:dpk
    !
    IMPLICIT NONE
    !ARG!
    INTEGER                         :: nmx
    REAL(dpk),    DIMENSION(nmx,nmx):: mx
    CHARACTER(len=*)                :: smx
    CHARACTER(len=*)                :: sf
    !LOC!
    INTEGER                         :: i,j
    


    WRITE(*,*) '# print_mx::', smx
    WRITE(*,*) '..........................................'

    DO i = 1, nmx 
       IF( sf == 'f') THEN  
!          WRITE(*,'A1,(5F8.1),A2') "|",( mx(i,j), j =1,nmx),"|"

          WRITE(*,'(30F10.5)') ( mx(i,j), j =1,nmx)
       ELSE IF( sf =='e' ) then
!          WRITE(*,'A1,(5E8.1),A2') "|",( mx(i,j), j =1,nmx),"|"

          WRITE(*,'(30E30.5)') ( mx(i,j), j =1,nmx)

       ENDIF
    ENDDO
    WRITE(*,*) '..........................................'

    RETURN
  END SUBROUTINE print_mx
!!!!!!!!!!!!!!!!!! OPEN FILES IN DATA MODE

  SUBROUTINE DATAFILE(N, L, FILENAME, TYPE)

    IMPLICIT NONE
    !ARG!
    INTEGER                                :: N
    INTEGER                                :: L
    CHARACTER(LEN=*)                       :: FILENAME
    CHARACTER(LEN=6), INTENT(in), OPTIONAL :: TYPE
    !LOC!
    CHARACTER( LEN = 30  )                 :: DATA_FILE
    CHARACTER( LEN = 15  )                 :: SL
    !EXE!


    WRITE(SL,'(I15)') L


    IF(L.GE.0) THEN

       IF(PRESENT(TYPE)) THEN

          DATA_FILE = "dat/"//FILENAME//TRIM(ADJUSTL(SL))//"." &
                           //TRIM(ADJUSTL(TYPE))               &
                           //".dat"
       ELSE
          DATA_FILE = "dat/"//FILENAME//TRIM(ADJUSTL(SL))//".dat"
       ENDIF
    
    ELSE IF(l.LT.0) THEN

       DATA_FILE = "dat/"//FILENAME//".dat"
    ENDIF

    !

    OPEN(n, file=data_file,form='unformatted',access='sequential')
    
    RETURN

  END SUBROUTINE DATAFILE
  !
!
!
! ascii mode
!
  SUBROUTINE asciifile(n, l, filename, type)
    !
    IMPLICIT NONE
    !
    !ARG!
    INTEGER                                :: N
    INTEGER                                :: L
    CHARACTER(LEN=*)                       :: FILENAME
    CHARACTER(LEN=6), INTENT(in), OPTIONAL :: TYPE
    !
    CHARACTER( LEN = 30  )                 :: ASCII_FILE
    CHARACTER( LEN = 15  )                 :: SL
    !
    !........
    

    WRITE(SL,'(I15)') L


    IF(L.GE.0) THEN

       IF(PRESENT(TYPE)) THEN

          ASCII_FILE = "out/"//FILENAME//TRIM(ADJUSTL(SL))//"." &
                           //TRIM(ADJUSTL(TYPE))               &
                           //".out"
       ELSE
          ASCII_FILE = "out/"//FILENAME//TRIM(ADJUSTL(SL))//".out"
       ENDIF
    
    ELSE IF(l.LT.0) THEN

       ASCII_FILE = "out/"//FILENAME//".out"
    ENDIF

    OPEN(n,file=ascii_file)
    RETURN
  END SUBROUTINE ASCIIFILE
  !
  !
  !
  SUBROUTINE afile(n,l1, l2, filename, smode)
    !
    IMPLICIT NONE
    !ARG!
    INTEGER                :: n
    INTEGER                :: l1
    INTEGER                :: l2
    CHARACTER(LEN=*)       :: filename
    CHARACTER( LEN = 30  ) :: asciifile
    CHARACTER( LEN = 6  )  :: smode
    !LOC!
    CHARACTER( LEN = 15  )  :: sl1,sl2
    !EXE!
    !
    !........
    
    WRITE(sl1,'(I15)') l1
    WRITE(sl2,'(I15)') l2

    asciifile  = "out/"//FILENAME & 
                      //TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                       //"."//TRIM(ADJUSTL(smode))              &  
                      //".out"

    OPEN(n,file=asciifile)
    WRITE(*, *) "# utils::dfile                 file = ",filename 
    RETURN
  END SUBROUTINE AFILE
  !
  !
  !
  !
  SUBROUTINE tafile(N, L, FILENAME, TYPE)

    IMPLICIT NONE
    !ARG!
    INTEGER                                :: N
    INTEGER                                :: L
    CHARACTER(LEN=*)                       :: FILENAME
    CHARACTER(LEN=6), INTENT(in), OPTIONAL :: TYPE
    !LOC!
    CHARACTER( LEN = 30  )                 :: DATA_FILE
    CHARACTER( LEN = 15  )                 :: SL
    !EXE!


    WRITE(SL,'(I15)') L


    IF(L.GE.0) THEN

       IF(PRESENT(TYPE)) THEN

          DATA_FILE = "tdat/"//FILENAME//TRIM(ADJUSTL(SL))//"." &
                           //TRIM(ADJUSTL(TYPE))               &
                           //".dat"
       ELSE
          DATA_FILE = "tdat/"//FILENAME//TRIM(ADJUSTL(SL))//".dat"
       ENDIF
    
    ELSE IF(l.LT.0) THEN

       DATA_FILE = "tdat/"//FILENAME//".dat"
    ENDIF

    !

!    OPEN(n, file=data_file,form='unformatted',access='sequential')
    
    OPEN(n,file=data_file)

    RETURN

  END SUBROUTINE tafile
    
  !
  !
  SUBROUTINE dfile(n, l1, l2, filename, smode)
    !
    IMPLICIT NONE
    !ARG!
    INTEGER                :: n
    INTEGER                :: l1
    INTEGER                :: l2
    CHARACTER(LEN=*)       :: filename
    CHARACTER( LEN = 30  ) :: datafile
    CHARACTER( LEN = 6  )  :: smode
    !LOC!
    CHARACTER( LEN = 15  )  :: sl1,sl2
    !EXE!
    !
    !........
    
    WRITE(sl1,'(I15)') l1
    WRITE(sl2,'(I15)') l2

    datafile  = "dat/"//FILENAME & 
                      //TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                       //"."//TRIM(ADJUSTL(smode))              &  
                      //".dat"

    OPEN(n,file=datafile, form='unformatted',access='sequential')
    PRINT*,"# utils::dfile: file = ",datafile
    RETURN
  END SUBROUTINE DFILE
!
!
!
SUBROUTINE check(n,m,k) 
  integer n,m,k
  if(n/=m) then 
     write(*,*) "# check: Inconsistency problem : ", k
     write(*,*) "#                            m = ", n
     write(*,*) "#                            n = ", m
     stop
  endif
!
END SUBROUTINE check
!...........................



  ! Performs the MX*V multiplication 
  !
  !    (1)  only the upper part of the square symmetric matrice (dim x dim) is provided 
  !
  !    (2)  the Bsp matrices are diagonally banded with width 2*width+1
  
  ! The above are exploited below

SUBROUTINE banded_mat_v_mul(mx, v, width, type)  
    !ARG!
    REAL(DPK), DIMENSION(:,:), INTENT(in)    :: mx
    REAL(DPK), DIMENSION(:), INTENT(inout)   :: V
    INTEGER                                  :: width
    INTEGER,         INTENT(in), OPTIONAL    :: type
    REAL(DPK), ALLOCATABLE, DIMENSION(:)     :: C
    !LOC!
    INTEGER I, J
    INTEGER fac, ierror
    !EXE!
    if(present(type)) then
     fac = type
     else
        fac = 1
    endif
    ierror = 1
    !  check if dimensions conform for mx(dim_mx,dim_mx) * vr(dim_v)

    CALL check(SIZE(v),SIZE(mx,1), ierror)
    ALLOCATE(C(SIZE(v)))
    
    c = 0.0_dpk
    DO i = 1, SIZE(v)
       !
       c(i) = c(i) + mx(i,i) * v(i)  ! diagonal part
       !
       lower_off_diagonal:DO j = MAX(1, i- width + 1), i - 1 
          c(i) = c(i) + fac * mx(j,i) * v(j) 
       ENDDO lower_off_diagonal
       !
       upper_off_diagonal:DO j = i + 1, MIN(i + width - 1, SIZE(v))
          c(i) = c(i) + mx(i,j) * v(j) 
       ENDDO upper_off_diagonal
    ENDDO
     
    v = c

    deallocate(c)
    RETURN
    
  END SUBROUTINE banded_mat_v_mul


END MODULE utils

