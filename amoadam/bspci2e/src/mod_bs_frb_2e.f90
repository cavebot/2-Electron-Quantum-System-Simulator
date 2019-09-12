
!
! modules:     cnstr_data:
!                    sub:  set_space: set_overlap_file(1), set_config(2), set_cnstr_array(3)
!                ite_data:
!                    sub: set_ite_data   
!


MODULE filenames
  !
  IMPLICIT NONE
  !
  CHARACTER(len=30) :: file_target_energies = 'dat/en1e.dat'
  !
END MODULE filenames


MODULE bs_frb_2e

  USE PRECISION, only:dpk
  IMPLICIT NONE
  PUBLIC
  
  
CONTAINS 


SUBROUTINE READ_TARGET_ENERGIES(en)
  !mod!
  USE PRECISION,  ONLY : DPK
  USE filenames,   only: file_target_energies
  !
  IMPLICIT NONE
  !arg!
  REAL(DPK), DIMENSION(:,:), POINTER :: en
  !loc!
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nmin, nmax 
  INTEGER                            :: lmin, lmax
  INTEGER                            :: nfile
  INTEGER                            :: i, j
  INTEGER                            :: dim_1, dim_2
  !exe!
  WRITE(*,'(a60)') 'bs_frb_2e::read_target_energies  in'
  NFILE = 46

  ! read and set up dimensions for target-energy matrix

  OPEN(nfile,file=file_target_energies)
  
  READ(nfile,*) lmin, lmax

  dim_1 = lmax-lmin+1

  allocate(nmin(dim_1))
  allocate(nmax(dim_1))

  READ(nfile,*) nmax(1), nmin(1)
  CLOSE(nfile)

  dim_2 = nmax(1)-nmin(1) + 1
 
  WRITE(*,'(a60,i10)') 'dim_1  = ', dim_1
  WRITE(*,'(a60,i10)') 'dim_2  = ', dim_2

  ALLOCATE(en(dim_1, dim_2))
  en = 0.0_dpk

  ! read entries for the target-energy matrix  

  OPEN(nfile,file=file_target_energies)


  READ(nfile,*) LMIN, LMAX  
  DO I = 1, dim_1

     READ(NFILE,*) nmax(i), nmin(i)

     IF(dim_2.NE.(nmax(i)-nmin(i)+1)) THEN
        WRITE(*,*) '& read_target_energies::                       warning: ' 
        WRITE(*,'(a60,i10)') '& read_target_energies:: nof eigenstates differ for i = ', i 
        WRITE(*,'(a60,i10)') '& read_target_energies::                    dim_1(i)  = ', nmax(i)-nmin(i)+1
        WRITE(*,'(a60,i10)') '& read_target_energies::                    dim_1(1)  = ', dim_2
     ENDIF

     DO j = 1, nmax(i) - nmin(i) + 1 
        READ(NFILE,'(E25.14)') EN(I,J) 
     ENDDO
  ENDDO
  CLOSE(NFILE)
  WRITE(*,'(a60)') 'bs_frb_2e::read_target_energies  in'
  RETURN
END SUBROUTINE READ_TARGET_ENERGIES
!
!
!
!

REAL(DPK) FUNCTION bvalue(t, bcoef, n, k, x, jderiv)
  !
  IMPLICIT NONE
  !
  REAL(DPK),      DIMENSION(n+k) :: t     
  REAL(DPK),      DIMENSION(n)   :: bcoef
  INTEGER                        :: n
  INTEGER                        :: k
  REAL(DPK)                      :: x
  INTEGER                        :: jderiv
  !loc!
  REAL(DPK),        DIMENSION(n) :: aj(n)
  REAL(DPK),        DIMENSION(n) :: dl(n)
  REAL(DPK),        DIMENSION(n) :: dr(n)
  REAL(DPK)                      :: fkmj
  !
  INTEGER   mflag, i, km1
  INTEGER   j, jcmin, imk, jcmax, ilo, jj, nmi, jc, kmj
  EXTERNAL interv  
  !exe!
     !
  bvalue = 0.0_dpk

  last_knot:IF(x.EQ.t(n+k)) THEN      !P(R) = c_n
     bvalue = bcoef(n)
  END IF last_knot
      

  IF(jderiv .GE. k)                go to 99
  CALL interv (t,n+k,x,i,mflag)

  IF(mflag .NE. 0)                 go to 99

  km1 = k - 1
  IF(km1 .GT. 0)                  go to 1
  bvalue = bcoef(i)
  go to 99
     
 1    jcmin = 1
      imk = i - k
      IF(imk .GE. 0)                    go to 8
      jcmin = 1 - imk

      DO 5 j = 1,i
         dl(j) = x - t(i+1-j)
5        CONTINUE
         DO 6 j = i,km1
         aj(k-j)  = 0.0D0
         dl(j) = dl(i)
6        CONTINUE
                                      go to 10

8                                     DO 9 j = 1,km1
         dl(j) = x - t(i+1-j)
9        CONTINUE

10       jcmax = k
      nmi = n - i
      IF(nmi .GE. 0)                  go to 18
      jcmax = k + nmi
      DO 15 j = 1,jcmax
         dr(j) = t(i+j) - x
15       CONTINUE
         DO 16 j = jcmax,km1
          aj(j+1) = 0.0D0
          dr(j) = dr(jcmax)
16        CONTINUE
          go to 20

 18   do 19 j = 1,km1
             dr(j) = t(i+j) - x
19           CONTINUE

 20   do 21 jc = jcmin,jcmax
                aj(jc) = bcoef(imk + jc)
21              CONTINUE
                
      if(jderiv .eq. 0)             go to 30
      DO 23 j = 1,jderiv
         kmj = k - j
         fkmj = DFLOAT(kmj)
         ilo = kmj
         DO 22 jj = 1,kmj
            IF(dl(ilo)+dr(jj).EQ.0.0d0) &
                 &WRITE(*,*) 'hey', dl(ilo+1)+dr(jj-1), fkmj 
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
22          CONTINUE
23          CONTINUE

30          IF(jderiv .EQ. km1)           go to 39
            DO 33 j = jderiv+1,km1
               kmj = k - j
               ilo = kmj
               DO 32 jj = 1,kmj
                  aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
                  IF(dl(ilo)+dr(jj).EQ.0.0d0) &
                       & WRITE(*,*) ' HEY 1 ', dl(ilo+1)+dr(jj-1)
                  ilo = ilo - 1
32                CONTINUE

33                CONTINUE

39                bvalue = aj(1)

99                CONTINUE

   END FUNCTION



END MODULE bs_frb_2e
!
!
!
MODULE read_data

!  USE PRECISION, ONLY:dpk

  PUBLIC write_en, write_v_mx, read_v_mx, write_wf1e, read_v 


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! store matrix elements ( upper - diagonal )
!
  SUBROUTINE write_mx(l, mx, filename, mx_type)
    !
    USE PRECISION,  ONLY : dpk
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
    USE PRECISION,  ONLY : dpk
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

    USE PRECISION,  ONLY : dpk
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

  use precision,  only : DPK
  use UTILS,      only : ASCIIFILE
  !
  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:) :: V
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I
  integer dim 
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

  use precision,  only : DPK
  use UTILS,      only : ASCIIFILE
  !
  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:) :: V
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I, J
  integer dim 
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
SUBROUTINE write_v_mx(l, v, mx, filename)
  !mod!
  USE PRECISION,  ONLY: dpk
  USE UTILS,      ONLY: datafile, check
  !
  implicit none
  !arg!
  integer L
  real(DPK), dimension(:)   :: v
  real(DPK), dimension(:,:) :: mx
  character(LEN=*)          :: filename
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
  store_v_mx:DO I = 1, dim_1
     WRITE(ifile)  v(i)  
     WRITE(ifile)( mx(i,j), j = 1, dim_2 )
  ENDDO store_v_mx
  CLOSE(ifile)

  WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
  WRITE(*,'(a60)') " subroutine write_v_mx out.", filename
  RETURN
  
END SUBROUTINE write_v_mx
!.......................
SUBROUTINE read_v_mx(l, v, mx, filename)
  !mod!
  USE PRECISION,  ONLY : dpk
  USE UTILS,      ONLY : datafile
  !
  IMPLICIT NONE
  !arg!
  INTEGER L
  REAL(DPK), DIMENSION(:),   POINTER :: v
  REAL(DPK), DIMENSION(:,:), POINTER :: mx
  CHARACTER(LEN=*)                   :: filename
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
  get_v_mx:DO i = 1, dim_1
     READ(ifile)  v(i)  
     READ(ifile)( mx(i,j), j = 1, dim_2 )
  ENDDO get_v_mx
  CLOSE(ifile)

  WRITE(*, '(a2,1X,a58,i10)') "#","    l = ", l
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_1 = ", dim_1
  WRITE(*, '(a2,1X,a58,i10)') "#","dim_2 = ", dim_2
  WRITE(*,'(a60)') " subroutine read_v_mx out."
  RETURN
END SUBROUTINE read_v_mx


END MODULE READ_DATA

!eof
