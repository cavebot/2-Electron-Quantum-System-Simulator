!
!
!
PROGRAM dmxab
  !mod!
  USE PRECISION, only:dpk
  USE dmx_1e_matrix
  USE param
  USE set_grid
  USE ioroutines
  USE utils
  !loc!
  IMPLICIT NONE
  INTEGER n1, k1, lang, l, ndim, i, j, nm, MODE
  INTEGER NOUT, N1E, NRBB, NRWFB
  REAL(dpk) dr0, rmax1, energy
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: coef
  CHARACTER(LEN=3) GAUGE
  CHARACTER(LEN=100) ARGV
!///////////

  CALL input
  CALL r_grid

  nm = 5

  CALL GETARG(1, ARGV)   
  READ(ARGV,*) GAUGE


  MODE = 1
  IF(gauge.EQ.'v')     MODE = 0

  WRITE(*,*) '# hdmx1e::    gauge = ', gauge
  WRITE(*,*) '# hdmx1e::     mode = ', mode


  N1E   = 1
  NRWFB = 3
  NRBB  = 4
  NOUT  = 16


  OPEN(NOUT, FILE='out/hdmx1e.out')

!///

  ndim = MAXVAL(n2)

  ALLOCATE(coef(ndim,ndim))

  lmax = 2
      
  DO l = Lmin, Lmax - 1

     WRITE(*,*) '#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

     WRITE(*,*) '# read fxd radial wfs for l = ', l-1

         coef = 0.0D+00

         CALL D1EFILE(N1E, l-1) 

         READ(N1E) dr0, rmax1, n1, k1, lang

         call check_real_input( dr0,   rs(l) )
         call check_real_input( rmax1, rmax  )
         call check_int_input (  n1,   n2(l) )
         call check_int_input (  k1,   k     )
         call check_int_input (  lang, l-1   )

         ! read  E,C_E(i) 

         DO i = 1, n1 - 2
            READ(n1e) energy, ( coef(j,n1-i-1), j = 2, n1 - 1)
         END DO         
         CLOSE(N1E)


         !
         ! Calculation of T(i, la ; j,lb) = < B_i(r),la| P_{n_a,l_a} >
         ! i = 1,...,nb
         ! n_a = 1,2,..., nb-2
         !
         
         CALL rbbfile(nrbb, l-1) 
         CALL overlap(nrbb, l-1, coef, ndim)        ! store as T_G(j,n_a)
         CLOSE(nrbb)
         
         
         !Calculation of <P_{n_al_a}| T_G | B_j, lb> ,   lb = la + 1
         ! P_{n_a,l_a} fxd orbital (na = 1,..,nb-2)
         ! j = 1,...,nb

         CALL RWFBFILE(nrwfb, l-1, l, mode )         ! store as T_G(j,n_a)
         CALL dmx( nrwfb, l-1, l, coef, ndim, mode)
         CLOSE(nrwfb)


         !

         WRITE(*,*) '# read fxd radial wfs for l = ', l
         
         CALL D1EFILE(N1E, l) 

         READ(N1E) dr0, rmax1, n1, k1, lang

         call check_real_input( dr0,   rs(l+1) )
         call check_real_input( rmax1, rmax    )
         call check_int_input (  n1,   n2(l+1) )
         call check_int_input (  k1,   k       )
         call check_int_input (  lang, l       )
          
         ! read en, C(en)

         DO i = 1, n1 - 2
            READ(N1E) energy,( coef(j, n1 - i - 1), j = 2, n1 - 1 )
         END DO
         CLOSE(N1E)

!
!           < P(l) | D | B(l-1) >          
!

         ! Calculation of <P_{n_al_a}| T_G | B_j, lb> ,   la = lb + 1
         ! P_{n_a,l_a} fxd orbital (na = 1,..,nb-2)
         ! j = 1,...,nb


         CALL RWFBFILE(nrwfb, l, l-1, mode) 
         CALL dmx(nrwfb, l, l-1, coef, ndim, mode )

         CLOSE(nrwfb)

      END DO


      ! calculate overlap for l = Lmax

      l = Lmax

      coef = 0.0D+00
      CALL D1EFILE(N1E, l-1) 

      READ(N1E) dr0, rmax1, n1, k1, lang
      CALL check_real_input( dr0,   rs(l) )
      CALL check_real_input( rmax1, rmax  )
      CALL check_int_input (  n1,   n2(l) )
      CALL check_int_input (  k1,   k     )
      CALL check_int_input (  lang, l-1   )
      DO i = 1, n1 - 2
         READ(N1E) energy, (coef(j,n1-i-1),j=2,n1-1)          ! read E,C_E(i)
      END DO
      CLOSE(N1E)

         !
         ! Calculation of T(i, la ; j,lb) = < B_i(r),la| P_{n_a,l_a} >
         ! i = 1,...,nb
         ! n_a = 1,2,..., nb-2
         !

         CALL RBBFILE(NRBB, l - 1)
         CALL overlap(NRBB, l-1, coef, ndim)
         CLOSE(NRBB)

       END PROGRAM dmxab
       !eof
             
