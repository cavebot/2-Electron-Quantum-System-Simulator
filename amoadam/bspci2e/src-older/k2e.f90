!
!
!
PROGRAM k2ef !(x kmtx)

  USE PRECISION, only: dpk
  USE PARAM
  USE SET_GRID
  USE RWF
  USE NORM
  USE IOROUTINES
  USE GET_SPACE
  USE bs_frb_2e, ONLY: read_target_energies
  !
  IMPLICIT NONE
  CHARACTER*16 infile, config, enfile
  INTEGER i, j, m, n, ne, nst, ic, ie, npp, nnp, ncount, inpl
  INTEGER nd, ndmax
  INTEGER ls, lo, ntot, ncsmx, nconfig, nmx, info, nbsp
  INTEGER ntmp
  INTEGER l_min, l_max, lp, nvalence
  INTEGER index, index1
  INTEGER NOUT, NCFG, N1E, NH2E, NAFILE, NKFILE, NZFILE

  INTEGER, ALLOCATABLE,DIMENSION(:)      :: ipiv
  REAL(DPK) energy, zin, e1, alpha, vp
  REAL(DPK) rcond
  REAL(DPK) tmp
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: coe_bsp!, ehf
  REAL(DPK), DIMENSION(:,:),  POINTER    :: ehf
  REAL(dpk), DIMENSION(2)                :: det

!!!new implementation
  REAL(dpk),    ALLOCATABLE, DIMENSION(:,:) :: aa, bb, kkx   ! k = a^-1*b
  REAL(dpk),    ALLOCATABLE, DIMENSION(:,:) :: wb            ! wb= P_{ij}(R)
  REAL(dpk),    ALLOCATABLE, DIMENSION(:,:) :: ux, invAA
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:,:) :: sx, s1x

  REAL(dpk),    ALLOCATABLE,DIMENSION(:)    :: eth, zeta, phi
  INTEGER,      ALLOCATABLE,DIMENSION(:)    :: n1, l1,l2
  REAL(dpk),    ALLOCATABLE, DIMENSION(:)   :: work
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:)   :: zwork
  COMPLEX(dpk)                              :: uim
  CHARACTER(LEN=100)  ARGV

!!! READ THE ORBITAL ANGULAR MOMENTUM
!!!..............................


      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL         
      !
      OPEN(8, file='inp/k2e.inp')     ! to be removed (read by para == no)
      READ(8,*) npp                      !nof points
      CLOSE(8)
      !
      alpha = 0.0_dpk
       vp   = 0.0_dpk
       zin  = 1.0_dpk
       !...!

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
    CALL getnumberofchannels(ncfg, lo, ncsmx, nbsp)  ! get nof channels
    CALL set_space(ncsmx, nbsp)                   !allocate arrays
    

    WRITE(*, '(a60,i10)') 'nof b-splines nb  = ',  nbsp
    WRITE(*, '(a60,i10)') 'nof channels) ncs = ',  ncsmx
        
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

    !xx  read     1e- energies for core states (thresholds)  ehf(l,i)

    CALL read_target_energies(ehf)                 ! read 1e energies (in Ryd)

    l_max = SIZE(ehf,dim=1) - 1
    !    WRITE(*,'(<l_max+1>E20.12)') (ehf(lp,1), lp=1,l_max+1)
    WRITE(*,'(10E20.12)') (ehf(lp,1), lp=1,l_max+1)
    WRITE(*,'(a60,i10)') '           nof target    waves = ', SIZE(ehf,dim=1)
    WRITE(*,'(a60,i10)') '           nof target energies = ', SIZE(ehf,dim=2)


    !determine maximum number of continuum channels
  
    ndmax = 0
    DO i = 1, ncsmx
       IF(idcs(i).EQ.1) THEN
          ndmax = ndmax + 1
       END IF
    END DO
    
    nconfig = 0
    nconfig = SUM(ndi) ! = ntot
    
    
    WRITE(*, *) '# k2e::   maximum nof open channes    ndmax = ',  ndmax, ntot
    WRITE(*, *) '# k2e::   total   nof cfgs          nconfig = ',  nconfig, ncsmx
    !note that ntot = nconfig
    

    !......     expand energy array   ehf(l,j) ----> eth(i)    


  ALLOCATE( eth(ndmax) )             !core energies
  DO i = 1, ndmax
     eth(i) = ehf( lhf(i) + 1, nhf(i))
  END DO
  



!  WRITE(*, *) ( eth(i), i = 1, ndmax)


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
     ALLOCATE(coe_bsp(nconfig,nd))
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

     ALLOCATE(    wb(nd, nd)     )           ! real quantities

     !
     ! read B-splines coefficients for each core channel
     !
     ! psi = S_{ij} C_ij * \phi_i(1) B_j(2) 


     READ(nh2e) ( (coe_bsp(n,m), n = 1, nconfig), m = 1, nd)

     
     WRITE(*, *) '#######################################################'
     WRITE(*,'(a60,i10)')    'energy index ie = ', ie
     WRITE(*,'(a60,E20.12)')               'E = ', energy
     WRITE(*,'(a60,i10)')  'nof open channels = ',  nd


     aa = 0.0_dpk         ! A - matrix        R(i,j)  --> A(i,j) * sin(w_i) + B(i,j) * cos(w_i)
     bb = 0.0_dpk         ! B - matrix           
     ux = 0.0_dpk         ! diagonal       A^-1 * R(i,j) -> sin(w_i) + K_ij cos(w_i)
     wb = 0.0_dpk         ! amplitudes matrix at the boundaries

     !     make ux diagonal
     DO i = 1, nd
        ux(i,i) = 1.0_dpk
     END DO

     index = 0 
     channels_in_final_state: DO ic = 1, nd

        nst    = 0
        ncount = 0
        coe    = 0.0_dpk

        WRITE(*,'(a60,i10)') '# solution index  ic = ',  ic

        all_cfg_series: DO i = 1, ncsmx                       ! all channels

           get_nof_open_channels: IF(idcs(i) == 1) THEN      ! continuum part of wf2e
              !ncount : open channel index out of ncsmx channel series
              IF( energy > eth(i) ) THEN
                 ncount = ncount + 1
                 WRITE(*,'(a19,1X,a21,I3,a5,1X,E15.7)') '# energy threshold', ' eth(',i,') = ', eth(i)
              ELSE 
                 CYCLE                
              ENDIF
           ENDIF get_nof_open_channels

           WRITE(*,'(a60,i10)') '# open channel index = ',  ncount
          
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

                                           !-----> i-th channel component [ehf(i),lhf(i),ll(i)] 
          coe( 1            ) = 0.0_dpk    !--------> zero b.c at the axis origin
          coe( 2:ndi(i) + 1 ) = coe_bsp( nst + 1:nst + ndi(i), ic)
                                                               !  
                                                               !_> ic-th solution
          nst = nst + ndi(i)

          index1 = 0
          normalize_open_channels:IF(idcs(i) == 1) THEN 

             index1 = index1 + 1
             index  = index  + 1
             
             CALL cal_wf( ll(i), coe, npp)          !compute 'outer' radial 1e-wf [ll(i)]  

             WRITE(*,'(a60,i10)') '# channel component index i = ',  i
             WRITE(*,*) ie, ic, i, ndi(i) 
             WRITE(*,*) nhf(i),lhf(i),ll(i), eth(i)
             WRITE(*,*) coe(ndi(i)+1 ), wf(npp)

!            write(*,*) ( wf(nnp), nnp = 1, npp, 20)
             
             IF(ncount.LE.nd) THEN         ! Normalize outer solutions
                ! 
                n1(ncount) = nhf(i)
                l1(ncount) = lhf(i)
                l2(ncount) =  ll(i)
                e1         = eth(i)  ! core threshold energy

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


                !ic        --> i-th solution
                !ncount    --> j-th channel 

                wb(ic,ncount)  =   coe(ndi(i)+1)                 
                aa(ic,ncount)  =   an * dcos(w1)             !  a(i,j)
                bb(ic,ncount)  = - an * dsin(w1)             !  b(i,j) 
                zeta(ncount)   =   aa2                       !  z(rmax)
                phi(ncount)    =   ph2                       !  w(rmax)

             END IF

          END IF normalize_open_channels
       END DO all_cfg_series
    END DO channels_in_final_state


    DEALLOCATE(coe_bsp)

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


    WRITE(*,'(a60,E25.14)') 'check inverse of A: A/A = ', MATMUL(invAA, aa)



    !,,,, A file (A^-1 stored)

   WRITE(NAFILE, '(1e15.7,2x,i5)') energy, nd
   DO i = 1, nd
      WRITE(NAFILE, '(4i5)') i, n1(i),l1(i),l2(i)
   ENDDO
   DO ic = 1, nd
      WRITE(NAFILE, '(5e15.7)')   ( aa(i,ic), i = 1, nd)     ! A^-1 - matrix in AMX-L.DAT
   END DO

   !,,,,,, Z  file  (zeta,phi stored)

   WRITE(NZFILE, '(1e15.7,2x,i5)')  energy, nd
   DO i = 1, nd
      WRITE(NZFILE, '(4i5,2E25.14)') i, n1(i),l1(i),l2(i), e1, energy-e1
   ENDDO
   DO ic = 1, nd
      WRITE(NZFILE, '(5e15.7)')   ( wb(i,ic), i = 1, nd)     ! amplitude - matrix in wb-L.DAT
   END DO
   !WRITE(NZFILE, '(5e15.7)')      ( zeta(i), i = 1, nd)   ! zeta(R) for each channel
   !WRITE(NZFILE, '(5e15.7)')      (  phi(i), i = 1, nd)   !  phi(R) for each channel


   kkx = MATMUL(aa, bb)         ! calculate K-matrix      u = pi * K = a^(-1) * b


   !,,,,, K file (K matrix stored) 

   WRITE(NKFILE, '(1e15.7,2x,i5)') energy, nd
   DO ic = 1, nd
      WRITE(NKFILE, '(5e15.7)') ( kkx(i,ic), i = 1, nd)   ! pi*K matrix in KMX-L.DAT
   END DO


   DEALLOCATE(bb)


   ! S-matrix calculation


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

!   WRITE(*,'(a60,E20.14)') 'check inverse of C=1+iK : C/C = ',  MATMUL( ( ux + uim * kkx ), s1x)


   DEALLOCATE( zwork )
   DEALLOCATE( ipiv  )

!!! calculate   S-matrix S = ( 1 - i * u ) / (1 + i * u )     ( u = pi * K ) 


   sx = MATMUL( ( ux - uim * kkx ), s1x)


   WRITE(NOUT,*) '# k2e::  en(ie) = ', ie, energy

   s_matrix_squared: DO i = 1, nd
      rcond = SUM( ABS(sx(i,:))**2 )
      WRITE(*,*) '# k2e::  |S(i)|^2 = ', i, rcond
  END DO s_matrix_squared

     DEALLOCATE(n1,l1,l2)
     DEALLOCATE(zeta)
     DEALLOCATE(phi)
     DEALLOCATE( wb)           ! real quantities
     DEALLOCATE( aa)           ! real quantities
     DEALLOCATE( invAA)
     DEALLOCATE( kkx )
     DEALLOCATE( ux )
     DEALLOCATE( sx )          ! complex quantities  
     DEALLOCATE( s1x ) 
  END DO final_2e_states
  
  CLOSE(NOUT)

END PROGRAM k2ef
!eof

