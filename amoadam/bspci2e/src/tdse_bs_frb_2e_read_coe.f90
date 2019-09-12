
PROGRAM coef_read

  USE PRECISION
  USE IO
  IMPLICIT NONE
  REAL(dpk),   ALLOCATABLE, DIMENSION(:,:) :: en 
  REAL(dpk),   ALLOCATABLE, DIMENSION(:)   :: y
  INTEGER,     ALLOCATABLE, DIMENSION(:  ) :: ne_L
  REAL(dpk) i0, omeg, tau, ecut, t, tim
  INTEGER nmax, n_states, ie
  INTEGER Lmax, nL, L
  INTEGER I,J, nt, ncoe
  CHARACTER(LEN=25) coe_t
  
  !... executable statements

  ncoe = 1

  !
  nmax = 1787

  OPEN(ncoe,file="tdat/coe.dat",form='unformatted',access='sequential')

  READ(ncoe) i0, omeg, tau, nL, ecut, tim
  Lmax = nL-1
  ALLOCATE( ne_L(0:Lmax) )
  ALLOCATE(   en(nmax, 0:Lmax) )
  
  READ(ncoe) ( ne_L(l), l = 0, lmax), n_states
  DO l = 0, lmax
     READ(ncoe) ( en(j,l), j = 1, ne_L(l))        !relative to he+(1s)
  END DO
  ALLOCATE(y(2*n_states))
  READ(ncoe) (y(i), i = 1, 2*n_states)         ! save 
  CLOSE(ncoe)
  
!end_read
  WRITE(*,*) i0, omeg, tau, nL, ecut, tim
  WRITE(*,*) lmax, n_states
  
  WRITE(*,*) "full vector"
  DO i = 1, n_states     
     WRITE(*, *)  i, i+n_states, y(i), y(i+n_states) 
  ENDDO

  nt = 0
  DO L = 0, lmax
     WRITE(*,*) L, ne_L(L)
     DO ie = 1, ne_L(L)        
        WRITE(*,'(2I6,1X,1E15.6,1x,1F12.4)') L, ie, y(ie + nt)**2+y(ie+nt+n_states)**2, en(ie,L) 
     ENDDO
      nt = nt + ne_L(L)
   ENDDO

   DEALLOCATE(y,ne_L,en)

  
END PROGRAM coef_read
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!EOF


