!
!
!  T_G(i, la ; j,lb) = < B_i(r),la| T_G | B_j, l >  '
!
!    i = 1,..., nb
!    j = 1,..., nb
!
!    la = 0,..,lmax
!    lb = 0,..,lmax
!
!    G = L (length), V (velocity)
! 

PROGRAM hdmx1eb

  USE PRECISION, only:dpk
  USE dmx_1e_matrix
  USE param
  USE set_grid
  USE ioroutines
  USE utils

  IMPLICIT NONE
  INTEGER n1, k1, lang, l, ndim, i, j, lf, nm, mode
  INTEGER NOUT, NRWFBB
  REAL(dpk) dr0, rmax1, energy
  CHARACTER(LEN=100) ARGV,GAUGE

!!!...............................


  CALL GETARG(1, ARGV)   
  READ(ARGV,*) GAUGE

  WRITE(*,*) '# hdmx1e::      1e dmx in gauge = ', gauge

  MODE = 1
  IF(gauge.EQ.'v')     MODE = 0     


  NRWFBB = 3
  NOUT  = 16


  call input
  call r_grid

  nm = 5

  OPEN(NOUT, FILE='out/hdmx1eb.out')

  ndim = MAXVAL( n2 )
  
  WRITE(*,*) '# dmx1eb::                      Lmin = ', Lmin
  WRITE(*,*) '# dmx1eb::                      Lmax = ', Lmax

!!!......................  

  DO l = Lmin, Lmax - 1
     WRITE(*,*) '# dmx1eb::                     L = ', l

     CALL RWFBBFILE(NRWFBB, l-1, mode )
     CALL     dmxbb(NRWFBB, l-1, l, ndim, mode)     
     CLOSE(NRWFBB)
         
  END DO

END PROGRAM hdmx1eb
!EOF
