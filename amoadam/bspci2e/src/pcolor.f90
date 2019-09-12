!
!
! transform data to matlab pcolor input
!
! x(n)   -> vector
! y(m)   -> vector
! p(n,m) -> matrix
!
!
PROGRAM write_pcolor_format
  
  USE PRECISION
  USE io
  IMPLICIT NONE
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: kx,ky,kz
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: p
  INTEGER nmax
  INTEGER i,j
  !
  
  
  ! read data

  nmax = 121380
  ALLOCATE(kx(nmax),ky(nmax),kz(nmax),p(nmax,nmax))

  kx = 0.0_dpk
  ky = 0.0_dpk
  kz = 0.0_dpk
  p  = 0.0_dpk

!xy
 
  OPEN(1,file='pad/pcolor.dat')
  OPEN(2,file='pad/kx.dat')
  OPEN(3,file='pad/ky.dat')
  OPEN(4,file='pad/p.dat')
  DO i = 1, nmax
     READ(1, '(4E20.8)') kx(i), ky(i), kz(i)
     WRITE(2,*) kx(i)
     WRITE(3,*) ky(i)

     WRITE(4,'(1188E20.8)') (p(i,j), j = 1, nmax)
  ENDDO
  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)

  DEALLOCATE(kx,ky,kz,p)

!!%!  READ(*,*) nmax
!!%!  WRITE(*,*)"# nmax = ", nmax 
!!%  nmax = 288
!!%  ALLOCATE(kx_xy(nmax),ky_xy(nmax),pxy(nmax,nmax))
!!%
!!%  kx_xy = 0.0d+00
!!%  ky_xy = 0.0d+00
!!%  pxy  = 0.0d+00
!!%
!!%!xy
!!% 
!!%  OPEN(1,file='tmp/kxy.dat')
!!%  OPEN(2,file='tmp/kx-xy.dat')
!!%  OPEN(3,file='tmp/ky-xy.dat')
!!%  OPEN(4,file='tmp/pxy.dat')
!!%  DO i = 1, nmax
!!%     READ(1, '(3E20.8)') kx_xy(i), ky_xy(i), pxy(i,i)
!!%     WRITE(2,*) kx_xy(i)
!!%     WRITE(3,*) ky_xy(i)
!!%     WRITE(4,'(297E20.8)') (pxy(i,j), j = 1, nmax)
!!%  ENDDO
!!%  CLOSE(1)
!!%  CLOSE(2)
!!%  CLOSE(3)
!!%  CLOSE(4)
!!%
!!%  DEALLOCATE(kx_xy,ky_xy,pxy)
!!%
!!%  nmax = 88
!!%  ALLOCATE(kx_xz(nmax),kz_xz(nmax),pxz(nmax,nmax))
!!%  kx_xz = 0.0d+00
!!%  kz_xz = 0.0d+00
!!%  pxz   = 0.0+00
!!%
!!%!xz
!!%  OPEN(1,file='tmp/kxz.dat')
!!%  OPEN(2,file='tmp/kx-xz.dat')
!!%  OPEN(3,file='tmp/kz-xz.dat')
!!%  OPEN(4,file='tmp/pxz.dat')
!!%  DO i = 1, nmax
!!%     READ(1, '(3E20.8)') kx_xz(i), kz_xz(i), pxz(i,i)
!!%     WRITE(2,*) kx_xz(i)
!!%     WRITE(3,*) kz_xz(i)
!!%     WRITE(4,'(1026E20.8)') (pxz(i,j), j = 1, nmax)
!!%  ENDDO
!!%  CLOSE(1)
!!%  CLOSE(2)
!!%  CLOSE(3)
!!%  CLOSE(4)
!!%
!!%  DEALLOCATE(kx_xz,kz_xz,pxz)


END PROGRAM write_pcolor_format



