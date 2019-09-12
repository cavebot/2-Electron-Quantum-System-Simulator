!!!
MODULE utils

  PUBLIC asciifile, datafile, check, write_mx, print_mx, c3j
contains



! print out matrix
  SUBROUTINE print_mx(nmx, mx, smx, sf) 

    USE PRECISION, ONLY:dpk
    IMPLICIT NONE
    INTEGER nmx
    REAL(dpk), DIMENSION(nmx,nmx)::mx
    CHARACTER(len=*) smx
    CHARACTER(len=*) sf
    INTEGER i,j
    !
    !
    WRITE(*,*) '# print_mx::', smx
    WRITE(*,*) '..........................................'

    DO i = 1, nmx 
       IF( sf == 'f') THEN  
!          WRITE(*,'A1,(5F8.1),A2') "|",( mx(i,j), j =1,nmx),"|"

          WRITE(*,'(30F10.3)') ( mx(i,j), j =1, nmx)
       ELSE IF( sf =='e' ) then
!          WRITE(*,'A1,(5E8.1),A2') "|",( mx(i,j), j =1,nmx),"|"

          WRITE(*,'(30E10.2)') ( mx(i,j), j =1, nmx)

       ENDIF
    ENDDO
    WRITE(*,*) '..........................................'

    RETURN
  END SUBROUTINE print_mx
!!!!!!!!!!!!!!!!!! OPEN FILES IN DATA MODE

  SUBROUTINE DATAFILE(N, L, FILENAME)

    implicit none

    integer N, L
    character(LEN=*)       FILENAME
    character( LEN = 30  ) DATA_FILE
    character( LEN = 6  )  SL

!........
    
    write(SL,'(I6)') L
    DATA_FILE = "dat/"//FILENAME//trim(adjustl(SL))//".dat"

    open(n, file=data_file,form='unformatted',access='sequential')
  
    return

  END SUBROUTINE DATAFILE

!!!!!!!!!!!!!!!!!!!! OPEN FILES IN ASCII MODE
  SUBROUTINE ASCIIFILE(N, L, FILENAME)
    IMPLICIT NONE
    INTEGER N, L
    CHARACTER( LEN = * ) FILENAME
    CHARACTER( LEN = 30) ASCII_FILE
    CHARACTER( LEN = 6 ) SL
    !
    WRITE(SL,'(I6)') L
    ASCII_FILE = "dat/"//FILENAME//TRIM(ADJUSTL(SL))//".dat"
    !
    OPEN(n, file=ascii_file)

    WRITE(*, '(a2,1X,a50,1X,a20)') "#","file = ",ascii_file
    !
    RETURN
  END SUBROUTINE ASCIIFILE
!!!!!!!!!!!!!!!!!
!
  SUBROUTINE check(n,m,er) 

    INTEGER n,m,er

    IF(n/=m) THEN 
       WRITE(*,*) "#utils::check: Inconsistency problem : ", er
       WRITE(*,*) "#                                  m = ", n
       WRITE(*,*) "#                                  n = ", m
       STOP
    ENDIF
    !
  END SUBROUTINE check
  !
END MODULE utils

