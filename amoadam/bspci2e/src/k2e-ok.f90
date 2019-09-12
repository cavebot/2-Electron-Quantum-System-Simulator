! LOG NOTES    laan/iesl/2001 
!
!   **  THIS PROGRAM RENORMALIZE THE 2-E EIGENSTATES OBTAINED   **
!   **  BY THE H2E.F90 PROGRAM. NORMALIZATION IS DONE ACCORDING
!   **  THE K-MATRIX NORMALIZATION
!   **
!      WRITTEN BY : JIAN ZHANG           ( 1994 - 1996 ) 
!                   L. A.A. NIKOLOPOYLOS ( 2001 - 2002 )
!
!      JIAN ZHANG : 1994-1996 
!                   RUN ON IBM MACHINES 
!                   XLF90 COMPILER/ ESSL ROUTINES 
!
!      L. A.A. NIKOLOPOULOS : 2001-2002
!
!       SEPT/OCT/NOV2001 : STRONG MODIFICATION OF THE I/O FACILITIES
!                        : CHANGING COMPILER TO PGF90 
!                        : CHANGING LIBRARY  TO LAPACK/BLAS
!                        : CHANGING OS TO LINUX
!
!
!      16022002          : MERGING THE KMTX.INP, H2E.INP FILE IN ONE
!                          NAMED 'H2E.INP' 
!                        : CREATING A NEW OUTPUT FILE FOR WARNING 
!                          INFORMATION NAMED ' KMTX.WARN'
!  
!                      
!       INPUT  FILES     : H2E.INP, CFG-L.INP, dat/en1e.dat
!       OUTPUT FILES     : KMTX.OUT, KMTX.WARN
!                          AMXS.DAT, KMXS.DAT, ZMXS.DAT
!
!        TO BE DONE 
!        * Wed Sep 12 12:05:28 EEST 2001
!
!
!       (-1) get NCS parameter from CFG-L.INP file. use set_space module
!            nbs, nbsplines parameters should read by input files.
!
!       (0)  Modified to allow automatic selection of the
!            I/O files for the calculations of 
!            A-MATRIX, K-MATRIX.
!      
!       (1)  Modified to allow dynamic memory allocation
!       (2)  For cpu-time limitations, might be avoided the  
!            calculation of  S - matrix, at the end of program
!
!###########################################################################



PROGRAM k2e !(x kmtx)

  USE PRECISION, only: dpk
  USE PARAM
  USE SET_GRID
  USE RWF
  USE NORM
  USE IOROUTINES
  USE GET_SPACE

  IMPLICIT NONE
  CHARACTER*16 infile, config, enfile
  INTEGER i, j, m, n, ne, nst, ic, ie, npp, nnp, ncount, inpl
  INTEGER nd, ndmax
  INTEGER ls, lo, ntot, ncsmx, nconfig, nmx, info, nbsp
  INTEGER ntmp
  INTEGER LLMIN, LLMAX, LP, NVALENCE
  INTEGER index, index1
  INTEGER NOUT, NCFG, N1E, NH2E, NAFILE, NKFILE, NZFILE

  INTEGER, ALLOCATABLE,DIMENSION(:) :: ipiv
  REAL(DPK) energy, zin, e1, alpha, vp
  REAL(DPK) rcond
  REAL(DPK) tmp
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: er, ehf
  REAL(dpk), DIMENSION(2) :: det

!!!new implementation
  REAL(dpk),    ALLOCATABLE,DIMENSION(:)    :: eth, zeta, phi
  INTEGER,      ALLOCATABLE,DIMENSION(:)    :: n1, l1,l2
  REAL(dpk),    ALLOCATABLE, DIMENSION(:)   :: work
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:)   :: zwork
  REAL(dpk),    ALLOCATABLE, DIMENSION(:,:) :: aa, bb, kkx, ux, invAA
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:,:) :: sx, s1x
  COMPLEX(dpk) uim
  CHARACTER(LEN=100) ARGV

!!! READ THE ORBITAL ANGULAR MOMENTUM
!!!..............................


      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL         


      OPEN(8, file='inp/k2e.inp')     ! to be removed (read by para == no)
      READ(8,*) npp                      !nof points
      CLOSE(8)

      alpha = 0.0_dpk
       vp   = 0.0_dpk
       zin  = 1.0_dpk
!...

    NH2E       = 1
    NCFG       = 15
    NOUT       = 16
    N1E        = 18
    NAFILE     = 10            
    NKFILE     = NAFILE + 10
    NZFILE     = NAFILE + 20

!...


    OPEN( NOUT,     FILE= 'out/k2e.out')


    WRITE(*, *)  '# k2e::                  L =', INPL
    WRITE(*, *)  '# k2e::        EFFECTIVE Z = ', zin
    WRITE(*, *)  '# k2e::              alpha = ', alpha 
    WRITE(*, *)  '# k2e::                 vp = ', vp


    !           read configuration file   

  lo = inpl 

!  CALL get_nof_channels(ncfg, lo, ncsmx, nbsp)  ! get nof channels

  CALL GETNUMBEROFCHANNELS(NCFG, LO, NCSMX, NBSP)  ! get nof channels
  CALL set_space(NCSMX, NBSP)                      !allocate arrays


  WRITE(*, *) '# k2e::           nof b-slines basis    nb  = ',  nbsp
  WRITE(*, *) '# k2e::           nof series:(channels) ncs = ',  ncsmx
        
  ! get the channel series and related information
  ! nmax   ==   n bsplines
  ! ntot   ==   total number of 1-e states ---> matrix dimension

  CALL cxfin(ncfg,nhf,lhf,ll,nmin,nmax,noll,is,lo,ls,ncsmx,ntot, ndi,idcs)
  

!.......

  IF(inpl.NE.lo) THEN
     WRITE(NOUT,*) ' CFG-L.INP FILE : INCOSISTENT L. Exiting'
     STOP
  ENDIF

!.......

  ALLOCATE( ehf(max0(MAXVAL(lhf), MAXVAL(ll))+1, MAXVAL(nmax)) )
  

  !xx  read     1e- energies for core states (thresholds)  ehf(l,i)

3 FORMAT(6i5)
  
  OPEN(N1E, FILE='dat/en1e.dat')

  READ(N1E, *) LLMIN, LLMAX

  IF(LLMAX.GT.max0( MAXVAL(lhf), MAXVAL(ll)) + 1)  LLMAX = max0(MAXVAL(lhf), MAXVAL(ll)) + 1
  
 
  EHF = 0.0D+00
  DO LP = LLMIN, LLMAX
     READ(N1E, 3) NMX, NVALENCE
     DO IE = 1, NMX
        READ(N1E,'(E20.14)')   ehf(lp,ie)
     END DO
  END DO
  
  CLOSE(N1E)


 !determine maximum number of continuum channels
  
  ndmax = 0
  DO i = 1, ncsmx
     IF(idcs(i).EQ.1) THEN
        ndmax = ndmax + 1
     END IF
  END DO


  nconfig = 0
  nconfig = sum(ndi) ! = ntot


  WRITE(*, *) '# k2e::   maximum nof open channes       nd = ',  ndmax, ntot
  WRITE(*, *) '# k2e::   total   nof cfgs          nconfig = ',  nconfig, ncsmx
!note that ntot = nconfig


!......     expand energy array   ehf(l,j) ----> eth(i) 


  ALLOCATE( eth(ndmax) )
  DO i = 1, ndmax
     eth(i) = ehf( lhf(i) + 1, nhf(i))
  END DO
  

!  WRITE(*, *) ( eth(i), i = 1, ndmax)


!!!xxxxxxx     Take information of states of L angular momentum: energies, coefficients
!!!
!!!              'L.e.dat "       E+ < E < E++
!!!              'L.ee.dat'       E > E++

  CALL hfile(nh2e,"dat","h2e","bin",lo)

  READ(NH2E) ne            !  # number of energies included

  CALL input  
  CALL r_grid  

  WRITE(*, *) '# k2e::   nof 2e  states    ne =  ', ne

  ALLOCATE(wf(no))

  uim = (0.0D+00, 1.0_dpk)             !  uim**2 = - 1   : fantastic  i 


  CALL hfile(nafile,"dat","a2e","ascii",inpl)
  CALL hfile(nkfile,"dat","k2e","ascii",inpl)
  CALL hfile(nzfile,"dat","z2e","ascii",inpl)


  final_2e_states: DO ie = 1, ne

     READ(NH2E) energy, nd          !'nd' nof open channels at energy 'energy' (out of ndmax)
     ALLOCATE(    er(nconfig,nd) )
     ALLOCATE(  zeta(nd)         )
     ALLOCATE(   phi(nd)         )
     ALLOCATE(    aa(nd, nd)     )           ! real quantities
     ALLOCATE( invaa(nd, nd)     )
     ALLOCATE(    bb(nd, nd)     )
     ALLOCATE(   kkx(nd, nd)     )
     ALLOCATE(    ux(nd, nd)     )
     ALLOCATE(    sx(nd, nd)     )           ! complex quantities
     ALLOCATE(   s1x(nd, nd)     )
     ALLOCATE(    n1(nd)         )
     ALLOCATE(    l1(nd)         )
     ALLOCATE(    l2(nd)         )

! read B-splines coefficients for each core channel

     READ(NH2E) ( (er(n,m), n = 1, nconfig), m = 1, nd)

     
     WRITE(*, *) '############################################### 2-e state ie = ', ie
     WRITE(*,'(a14,1X,a26,I3,a5,1X,E15.7)') '# Total energy', ' E(',ie,') = ', energy
     WRITE(*,'(a49,1X,I5)') '#        nof open channels = ',  nd



     aa = 0.0_dpk         ! A - matrix        R(i,j)  --> A(i,j) * sin(w_i) + B(i,j) * cos(w_i)
     bb = 0.0_dpk         ! B - matrix           
     ux = 0.0_dpk
     
     DO i = 1, nd
        ux(i,i) = 1.0_dpk
     END DO

     index = 0 
     channels_in_final_state: DO ic = 1, nd

        nst = 0
        coe = 0.0_dpk
        ncount = 0

        WRITE(*,'(a28,1X,I5)') '# i-th open solution index = ',  ic
      
        all_cfg_series: DO i = 1, ncsmx                       ! all channels
           free_boundary_channels: IF(idcs(i) == 1) THEN      ! continuum part of wf2e

              !ncount : open channel index out of ncsmx channel series

              IF( energy > eth(i) ) THEN
                 ncount = ncount + 1
                 WRITE(*,'(a19,1X,a21,I3,a5,1X,E15.7)') '# energy threshold', ' eth(',i,') = ', eth(i)
              ELSE 
                 CYCLE                
              ENDIF
           ENDIF free_boundary_channels


          
!!!xxx     Calculate radial w.f. for outer e
!!!xxx  
!!          F_ij(r_2)     i-th solution, j-th channel component
!!         
!!          | F_ij > == Sum_j  c(i,j) * | f_j >     i = 1 ,...,  nd  i-th solution
!!                                                  j = 1 ,...., nd  j-th component
!!
!!          f_j(r_2)  :  uncorrelated  channels
!!          F_ij(r_2) :  correlated    channels
!!          c(i,j)    :  mixing coefficients       (1/r_12 - 1/r_2) 
!!




          coe( 1            ) = 0.0_dpk    !--------> zero b.c at the axis origin
                                            
                                           !-----> i-th channel component [ehf(i),lhf(i),ll(i)] 
                                           !
          coe( 2:ndi(i) + 1 ) = er( nst + 1:nst + ndi(i), ic)
                                                          !  
                                                          !_____> ic-th solution

          nst = nst + ndi(i)

          index1 = 0

          IF(idcs(i) == 1) THEN 

             index1 = index1 + 1
             index  = index  + 1
             
             CALL cal_wf( ll(i), coe, npp)          !compute 'outer' radial 1e-wf [ll(i)]  

!            write(*,*) ( wf(nnp), nnp = 1, npp, 20)
             
             IF(ncount.LE.nd) THEN         ! Normalize outer solutions

! 
                 n1(ncount) = nhf(i)
                 l1(ncount) = lhf(i)
                 l2(ncount) =  ll(i)


                e1 = eth(i)                        ! core threshold energy

                zin = 1.0_dpk                       ! effective charge
                IF(e1.GT.0.0D+00)   zin = 2.0_dpk

                WRITE(*,'(a49,1X,I5)') '# j-th component index = ',  ncount 
                WRITE(*,'(a18,1X,a22,I3,a5,1X,E15.7,1X,2I2 )')  &
                & '# inner 1-e:      ', '  e1(',ncount,') = ', e1, n1(ncount),l1(ncount)
                WRITE(*,'(a18,1X,a22,I3,a5,1X,E15.7,1X,I2)') &
                & '# outer 1-e:      ', '  e2(',ncount,') = ', energy-e1, l2(ncount)
                WRITE(*,'(a18,1X,a22,I3,a5,1X,G10.2)') '# effective charge', 'zeff(',ncount,') = ', zin

!
!compute  phase shift + normalization factor
!
!                F(i,j)  -->   A(i,j) * sin( w_i + w1(i) )
!                        -->   a(i,j) * sin( w_i ) + b(i,j) * cos(w_i)


                CALL nom( energy, ll(i), npp, r, wf, e1, zin, alpha, vp)

                aa(ic, ncount) =   an * dcos(w1)             ! A - matrix
                bb(ic, ncount) = - an * dsin(w1)             ! B - matrix 
                zeta(ncount)   =   aa2                       ! j(Rmax)
                phi(ncount)    =   ph2                       ! w_i(Rmax)


!!!                WRITE(*, *) ' NORMALIZATION CONSTANT AN = ', ABS(an)
!!!                WRITE(*, *) ' PHASE SHIFT            W1 = ', w1

             END IF
          END IF

       END DO all_cfg_series


    END DO channels_in_final_state


    DEALLOCATE(er)

    invaa = aa                    ! save aa


    ALLOCATE( work(64*nd) )
    ALLOCATE( ipiv(nd)  )

    CALL DGETRF( nd, nd, aa, nd, ipiv,  info)
    IF (INFO.EQ.0) THEN       !           Compute inverse of A       
       CALL DGETRI( nd, aa, nd, ipiv, work ,64*nd , info)
    ELSE
       WRITE(*,*) 'k2e:: matrix is singular. A^(-1) does not exist.'
       STOP
    ENDIF
    
    DEALLOCATE(work)



!   WRITE(*,*) 'A-1 :: '!WRITE(*,'(4e14.6)') ((AA(i,ic),i=1,nd),ic=1,nd) 


!!!    WRITE(NOUT,*) ' A * A^- 1 = ', MATMUL(invAA, aa)



!   WRITE(NOUT, *) ' A^-1  matrix in DAT/AMX-L.DAT'

   WRITE(NAFILE, '(1e15.7,2x,i5)') energy, nd
   DO i = 1, nd
      WRITE(NAFILE, '(4i5)') i, n1(i),l1(i),l2(i)
   ENDDO
   DO ic = 1, nd
      WRITE(NAFILE, '(5e15.7)')   ( aa(i,ic), i = 1, nd)     ! A^-1 - matrix in AMX-L.DAT
   END DO

   WRITE(NZFILE, '(1e15.7,2x,i5)')  energy, nd
   WRITE(NZFILE, '(5e15.7)')      ( zeta(i), i = 1, nd)   ! zeta(R) for each channel
   WRITE(NZFILE, '(5e15.7)')      (  phi(i), i = 1, nd)   !  phi(R) for each channel


   kkx = MATMUL(aa, bb)         ! calculate K-matrix      u = pi * K = a^(-1) * b


   DEALLOCATE(bb)



!!   WRITE(*, *) ' pi * K matrix in DAT/KMX-L.dat'

   WRITE(NKFILE, '(1e15.7,2x,i5)') energy, nd
   DO ic = 1, nd
      WRITE(NKFILE, '(5e15.7)') ( kkx(i,ic), i = 1, nd)   ! pi*K matrix in KMX-L.DAT
   END DO



   s1x = ux + uim * kkx


!s1x ---> s1x^{-1}


   ALLOCATE( zwork(64*nd) )

   CALL zgetrf(nd, nd, s1x, nd, ipiv, info)        ! Make the LU factorization

   IF(INFO == 0 ) THEN

   CALL zgetri(nd, s1x, nd, ipiv, zwork, 64*nd, info) 
   ELSE
      WRITE(*,*) 'k2e:: matrix is singular. S^(-1) does not exist.'
      STOP
   ENDIF

!!!   WRITE(*,*) ' S * S^(-1)  ',  MATMUL( ( ux + uim * kkx ), s1x)


   DEALLOCATE( zwork )
   DEALLOCATE( ipiv  )

!!! calculate   S-matrix S = ( 1 - i * u ) / (1 + i * u )     ( u = pi * K ) 


   sx = MATMUL( ( ux - uim * kkx ), s1x)


   WRITE(NOUT,*) '# k2e::  en(ie) = ', ie, energy

   s_matrix_squared: DO i = 1, nd
      rcond = SUM( ABS(sx(i,:))**2 )
      WRITE(NOUT,*) '# k2e::  |S(i)|^2 = ', i, rcond
  END DO s_matrix_squared


!!!      WRITE(NOUT,*) ((s1x(i,ic),i=1,nd),ic=1,nd) 
!!!      WRITE(NOUT,*) ((dreal(s1x(i,ic)),i=1,nd),ic=1,nd) 


!!!...........................

     DEALLOCATE(n1,l1,l2)
     DEALLOCATE(zeta)
     DEALLOCATE(phi)
     DEALLOCATE( aa)           ! real quantities
     DEALLOCATE( invAA)
     DEALLOCATE( kkx )
     DEALLOCATE( ux )
     DEALLOCATE( sx )          ! complex quantities  
     DEALLOCATE( s1x ) 
  END DO final_2e_states

  CLOSE(NOUT)

END PROGRAM k2e
!eof

