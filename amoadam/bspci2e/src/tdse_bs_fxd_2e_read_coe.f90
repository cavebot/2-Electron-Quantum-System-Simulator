
PROGRAM coef_read

  USE PRECISION
  USE IO
  IMPLICIT NONE
  REAL(dpk),   ALLOCATABLE, DIMENSION(:,:) :: en 
  REAL(dpk),   ALLOCATABLE, DIMENSION(:)   :: y
  INTEGER,     ALLOCATABLE, DIMENSION(:  ) :: ne_L
  REAL(dpk)                                :: time
  INTEGER                                  :: nmax, n_states, ntot
  REAL(dpk)                                :: i0, omeg, tau, ecut
  REAL(dpk)                                :: ry,iy,py
  INTEGER                                  :: Lmax, L
  INTEGER                                  :: i,j, nt, ncoe, ie

  
  !... executable statements

  ncoe = 1
  nmax = 2000

  !


  OPEN(ncoe,file="tdat/coe.dat",form='unformatted',access='sequential')

  !
  READ(ncoe) i0, omeg, tau, lmax, ecut


  ! get total nof states
  n_states = 0
  DO l = 0, lmax
     READ(ncoe) ie       
     n_states = n_states + ie
  END DO

  ! allocate
  ALLOCATE( ne_L(0:Lmax) )
  ALLOCATE( en(n_states, 0:Lmax) )
  ALLOCATE(  y(2*n_states) )



  REWIND(ncoe)
!!!!!!!!!!!!!!!!!!!!
  READ(ncoe) i0, omeg, tau, lmax, ecut
  DO l = 0, lmax
     READ(ncoe) ne_L(l), ( en(j,l), j = 1, ne_L(l))        !relative to he+(1s)
  END DO
  READ(ncoe) ntot 
  READ(ncoe) time, (y(i), i = 1, 2*n_states)         ! save 
!!!!!!!!!!!!!!!!!!!!



  IF(ntot.NE.n_states) THEN 
     PRINT*, "# something wrong in the input file, ntot not equal to sum of n_L, L = 0,1,..,Lmax"
     PRINT*, '#'  
     PRINT*, '#               ntot =', ntot  
     PRINT*, '#           nstates  =', n_states  
     PRINT*, '#   exit....'  
     STOP
  ENDIF


  CLOSE(ncoe)
  
!dump out

  WRITE(*,*) i0, omeg, tau, lmax, ecut, time
  WRITE(*,*) n_states
  
  WRITE(*,*) " vector: C_nl(t):  i, i + ntot, re_i = c(i), im_i = c(i+ntot):"
  DO ie = 1, n_states
     ry = y(ie)
     iy = y(ie+n_states)
     py = ry**2 + iy**2
     WRITE(*, *)  ie, ry,iy,py
  ENDDO
  WRITE(*,*) "&"
  WRITE(*,*) "&"
  WRITE(*,*) "density:  L, n, |C_nl(t)|^2 "
  WRITE(*,*) "&"
  nt = 0
  DO L = 0, lmax
     WRITE(*,*) L, ne_L(L)
     WRITE(*,*) '&   (nl,l) = (',ne_l(l),l, ')'
     DO ie = 1, ne_L(L)        
        WRITE(*,'(2I6,1X,1E15.6,1x,1F12.4)') L, ie, y(ie + nt)**2 + y(ie + nt + n_states )**2, en(ie,L) 
     ENDDO
      nt = nt + ne_L(L)
   ENDDO

   DEALLOCATE(y,ne_L,en)
  
END PROGRAM coef_read
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!EOF


