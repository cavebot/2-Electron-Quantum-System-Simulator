! LOG NOTES
!
!   **  THIS PROGRAM CALCULATES THE 2-E EIGENSTATES BY SOLVING   **
!   **  A LINEAR SYSTEM OF THE FORM A * X = B
!   **
!   **
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
!      16022002          : READ THE NUMBER OF CHANNELS FROM THE CONFIGURATIUON 
!                          FILE  CFG-L.INP 
!       INPUT FILES      : H2E.INP, CFG-L.INP
! 
!!!#####################################################################
MODULE kinds
  INTEGER, PARAMETER :: dpk = KIND(1.d0)
END MODULE kinds

MODULE cnstr_data
  USE kinds
  IMPLICIT NONE

  PUBLIC
  CHARACTER(len=16), ALLOCATABLE, DIMENSION(:) :: overlap
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ix, iy , ndi, ic, idcs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhf, lhf, ll, nmin, nmax, is, noll
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: cvec

  INTERFACE set_space
     MODULE PROCEDURE set_overlap_file, set_cnstr_array, set_config
  END INTERFACE

  PRIVATE set_overlap_file, set_cnstr_array, set_config
CONTAINS 

!!!##############################

  SUBROUTINE set_overlap_file(nl)

    IMPLICIT NONE
    INTEGER  nl

    ALLOCATE(overlap(nl))

    WRITE(*,*) '# h2e::set_overlap_file:: allocation for overlap matrices done. nl = ', nl


  END SUBROUTINE set_overlap_file

!!!###############################

  SUBROUTINE set_config(nd, nchl)

    IMPLICIT NONE
    INTEGER  nchl, nd

  
    ALLOCATE( ndi(nd)  )
    ALLOCATE( nhf(nd)  )
    ALLOCATE( lhf(nd)  )
    ALLOCATE( ll(nd)   )
    ALLOCATE( nmin(nd) )
    ALLOCATE( nmax(nd) )
    ALLOCATE( noll(nd) )
    ALLOCATE( is(nd)   ) 
    ALLOCATE( idcs(nd) )

    WRITE(*,*) '# h2e::set_config: set configurations done. nd = ', nd 
    
  END SUBROUTINE set_config

!!!###############################

  SUBROUTINE set_cnstr_array(nocnstr, ncsmx, nmx)
    IMPLICIT NONE
    INTEGER nocnstr, ncsmx, nmx

    ALLOCATE(ix(nocnstr)) 
    ALLOCATE(iy(nocnstr)) 
    ALLOCATE(ic(nocnstr)) 
    ALLOCATE(cvec(nmx, nocnstr))

    WRITE(*,*) '# h2e::set_cnstr_array: set constraints array done. noc = ', nocnstr 
  END SUBROUTINE set_cnstr_array

END MODULE cnstr_data

!!!#########################################

MODULE ite_data
  USE kinds
  IMPLICIT NONE
  PUBLIC
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: cx, bx
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: er, vec 
  REAL(DPK), ALLOCATABLE, DIMENSION(:) :: finit

  INTERFACE ite_space
     MODULE PROCEDURE set_ite_data
  END INTERFACE

  PRIVATE set_ite_data
CONTAINS

!!!#################################

  SUBROUTINE set_ite_data(ndim, nchl)
    IMPLICIT NONE
    INTEGER  ndim, nchl
!!!....................

    WRITE(*,*) '# h2e::set_ite:  allocation for energy matrix done. ', ndim, nchl 
    WRITE(*,*) '# h2e::set_ite:                              ndim = ', ndim 
    WRITE(*,*) '# h2e::set_ite:                              nchl = ', nchl 

    ALLOCATE(er(ndim,nchl))
    ALLOCATE(vec(ndim,nchl))
    ALLOCATE(finit(ndim))

  END SUBROUTINE set_ite_data


END MODULE ite_data

!!!###################################################
  PROGRAM H2E

    USE cnstr_data
    USE ite_data
    USE ioroutines

    IMPLICIT NONE
    
    INTEGER::  nm, j, i, n, m, nmx, icf, nx, ncol, nd, ne, nk, nkk,nrhs
    INTEGER::  ntot, ncsmx, lo, ls, no, nconfig, ntest, nchl, inpl
    INTEGER::  LP,IE, LLMIN,LLMAX, ncore, k , NDMAX
    INTEGER::  iflag, nocnstr, n0, lang, ndim, ns, nl, nvalence, NBSP
    INTEGER::  info, ITRY, nsymx
    INTEGER::  NOUT, NCFG, N1E, NR12, NHMXFILE, NBMXFILE
    INTEGER::  NH2E, NH2ESINGLE, NH2EDOUBLE
    REAL(DPK)::  ax, ay
    REAL(DPK)::  energy, dr0, rcond, e0, de, eMax, eTEST, EnAU 


    INTEGER, ALLOCATABLE, DIMENSION(:) :: ih, ipvt
    INTEGER :: lwork, ilaenv
    REAL(dpk), ALLOCATABLE, DIMENSION(:) :: lapack_work

    REAL(DPK), ALLOCATABLE, DIMENSION(:) :: aux
    REAL(DPK), ALLOCATABLE, DIMENSION(:) :: adiag, eth
    REAL(DPK), ALLOCATABLE, DIMENSION(:) :: fv1, fv2, ogp, er1
    REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: work, hx, ehf

    REAL(DPK), DIMENSION(2) :: det

    CHARACTER(LEN = 100) ARGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    NH2ESINGLE = 1                     !!! DO NOT CHANGE
    NH2EDOUBLE = 11
    NR12       = 3 
    NHMXFILE   = 7
    NBMXFILE   = 8
    NCFG       = 15
    NOUT       = 16 
    N1E        = 18

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL         
      CALL GETARG(2, ARGV)
      READ(ARGV,*) NE         

!!! READ THE ORBITAL ANGULAR MOMENTUM
!!!..............................

!!!..............................
    OPEN(NOUT, FILE='OUT/H2E.LOG')


    WRITE(*, *) '# h2e::                                 L =', INPL

!!!..............................


    OPEN(8, file='INP/H2E.INP')
    READ(8,*) e0, de                ! Read energy range parameters
    READ(8,*) nocnstr               ! Read degeneracy & no. constraints
    READ(8,*) ncore
    CLOSE(8)

!..............................


    eMax = e0 + ne * de 
        



    WRITE(*,*) '# h2e::        ei = ', e0,  ' a.u. = ', e0*27.211, ' eV' 
    WRITE(*,*) '# h2e::        de = ', de,  ' a.u. = ', de*27.211, ' eV'   
    WRITE(*,*) '# h2e::        ef = ', eMax,' a.u. = ', eMax*27.211, ' eV' 
    WRITE(*,*) '# h2e::        ne = ', ne    

    IF(INPL.LT.0) THEN 
       WRITE(*,*) '# h2e::              check'
       STOP
    ENDIF



    e0 = 2.0D+00 * e0            !go on with Ryd
    de = 2.0D+00 * de

!............


    !read configuration file


        lo = inpl

        ! get the number of chanell series from CFG-L.INP file

        CALL GETNUMBEROFCHANNELS(NCFG, LO, NDMAX, NBSP)

        ! allocate arrays for the quantum numbers of the channels 

        ! nchl is obsolete now. Is used for interface reasons,

        CALL set_space(NDMAX, NCHL)

        WRITE(*, *) '# h2e                          ndmax = ', NDMAX
        
        ! get the channel series and related information
        ! nmax   ==   n bsplines
        ! ncsmx  ==   ndmax
        ! ntot   ==   total number of 1-e states ---> matrix dimension

        CALL cxfin(ncfg, nhf,lhf,ll,nmin,nmax,noll,is,lo,ls, ncsmx, ntot, ndi, idcs)


!!!..............................

        IF(ls.EQ.1) THEN
           WRITE(*,*) '# h2e::                SINGLET SYMMETRY'
        ELSE IF(ls.EQ.3)   THEN
           WRITE(*,*) '# h2e                   TRIPLET SYMMETRY'
        ELSE
           WRITE(*,*) ' Allowed values of isl in r12.inp are 1 (Singlet), or 3 (triplet) '
        ENDIF
  
        IF(inpl.NE.lo) THEN
           WRITE(*,*) ' CFG-L.INP FILE : INCOSISTENT L. Exiting'
           STOP
        ENDIF

!!!..............................

        ns = MAXVAL(nmax)
        nl = max0(MAXVAL(lhf), MAXVAL(ll)) + 1

        ALLOCATE(ehf(nl, ns))


!!!...............................
!!!
!!! COMPUTE TOTAL NUMBER OF CONFIGURATIONS
!!! 
!!!   NDI(I)  :: NUMBER OF CONFIGURATIONS FOR CHANNEL I 

        nconfig = 0
        DO i = 1, ncsmx

           n = is(i) - 1

           DO j = 1, ndi(i)   

              n = n + 1
              nconfig = nconfig + 1

           END DO
        END DO


        ndim = MAX(ntot, nconfig + nocnstr)

!!!..................................


        WRITE(*, *) '# h2e::                      ncs  = ', nconfig
        WRITE(*, *) '# h2e::                   nconst  = ', nocnstr
        WRITE(*, *) '# h2e::                     ndim  = ', ndim


!!!....................................
!!!#     Read 1-e angular momenta and energies (en1e.dat)


3       FORMAT(6i5)

      OPEN(N1E, FILE='dat/en1e.dat')

      READ(N1E, *) LLMIN, LLMAX

      IF(LLMAX.GT.NL) LLMAX = NL

      DO  LP = LLMIN, LLMAX

         DO  IE = 1, NS

            EHF(LP,IE) = 0.0D+00
         ENDDO 

         READ(N1E, 3) NMX, NVALENCE
         DO  IE = 1, NMX
            READ(N1E,'(E20.14)') EHF(LP,IE)
         ENDDO
      ENDDO

      CLOSE(N1E)

!!!............................


      ALLOCATE( adiag(ncsmx) )
      ALLOCATE( hx(ndim, ndim) )
      ALLOCATE( bx(ndim, ndim) )
      ALLOCATE( eth(ncsmx) )

      WRITE(*,*) '# h2e::         , A, H, B allocation done. '

      hx    = 0.0D+00
      bx    = 0.0D+00
      adiag = 0.0D+00

      DO icf = 1, ncsmx

         SELECT CASE(idcs(icf))

         CASE(1)

            n = ll(icf) + 1
            
!!!  READ HAMILTONIAN MATRIX ON B-SPLINES  < B_I | H | B_J >
!!!                  AND
!!!  READ OVERLAP OF B-SPLINES  < B_I | B_J >


           CALL HMXFILE( NHMXFILE, ll(icf) ) 
           CALL BMXFILE( NBMXFILE, ll(icf) ) 

              
           READ(NHMXFILE) nmx
           READ(NBMXFILE) nmx

           DO i = 1, noll(icf)

              n   = is(icf) + i - 1
              nx  = is(icf) + noll(icf) - 1

              READ(NHMXFILE) ( hx(n,j), j = n, nx)
              READ(NBMXFILE) ( bx(n,j), j = n, nx)
               
               !                           build H for uncoupled channels 

              hx(n, n:nx) = hx(n,n:nx) + bx(n,n:nx) * ehf( lhf(icf) + 1, nhf(icf) )
              bx(n:nx, n) = bx(n,n:nx)

           END DO


           adiag(icf) = hx(nx, nx)

        CASE(0)

           DO i = 1, noll(icf)

              n       = is(icf) + i - 1
              nk      = i + nmin(icf) - 1
              hx(n,n) = ehf(lhf(icf) + 1, nhf(icf)) + ehf(ll(icf)+1,nk)
              bx(n,n) = 1.0D+00
              
           END DO

        END SELECT

                                !          write(*,*) 'adiag()=', icf, adiag(icf)

        CLOSE(NHMXFILE)
        CLOSE(NBMXFILE)

        eth(icf) = ehf( lhf(icf) + 1, nhf(icf))
                                !          write(*,*) 'eth=', icf, eth(icf)
     END DO


!!!..............................        



     DEALLOCATE(ehf)

!!!     WRITE(*,*) ' DIAGONAL ELEMENTS OF HAMILTONIAN'
!!!   !            write(*,'(3e12.5)') (hx(n,n),n=1,5)




!!!.........                      add V12 to the Hamiltonian        


!!
!!
!!    ndim < dimension of hxh matrix (calculated from hxh-mix.f) programm.
!!    or total number of configuration     ncf in    hxh.din should be more 
!!    than total number of configurations in dg2ec.din + no of constraints.
!!
!!       ncf(hxh.din) >= ncf(dg2ec.din) + ncostrnt
!!



     ALLOCATE(work(ndim,ndim)) 

     work = 0.0D+00

      WRITE(*,*) '# h2e::                    R12 allocation done.'
!!!....................



!!!.................... 
!!!     READ H_12 MATRIX ELEMENTS AND ADD TO FREE 
!!!     HAMILTONIAN
!!!     
!!!       H = H(1) + H(2) + H_12 
!!!


     CALL HR12FILE(NR12, LO)  


     DO i = 1, ntot
 
        READ(NR12) ( work(i,j),j = i, ntot)

        hx(i, :)         = hx(i, : )      +  work(i, : )
        hx(i:ntot, i )   = hx(i, i:ntot )

     END DO



     WRITE(*,*) '# h2e::                     r12 matrix (diagonal):'
     WRITE(*, '(3e12.5)')  ( work(n,n), n = 1, ntot, 41 )
     WRITE(*,*) '# h2e::          h(1)+h(2) matrix (diagonal part):'
     WRITE(*, '(3e12.5)')  ( hx(n,n),   n = 1, ntot, 41 )
        
     CLOSE(NR12)


!!..............................................
!!
!!
!!  WRITE(*,*)'WRITING HAMILTONIAN MATRIX :'
!!
!!..............................................

         WRITE(*,*) ' =========== HX MATRIX ============= '
      DO n = 1, ntot - 1,  ntot - 2

         WRITE(*,'(3e12.5)') ( hx(n,m), m = 1, ntot-1, ntot-2)
         
      END DO




!!!........................
        WRITE(*,*)  '# h2e::                           bx matrix:'
     DO n = ntot - 2, ntot

        WRITE(*,'(3e12.5)') (bx(n,m), m = n, ntot)
        
     END DO
   

!!!.......................

     
     ALLOCATE(ih(ndim))


     ih  = 0
     DO i = 1, ncsmx
        
        n = is(i) - 1
           
        DO j = 1, ndi(i)   

           n     =  n + 1
           ih(n) =  1

        END DO
     END DO



!     WRITE(*,*) '  NUMBER OF CONFIGURATIONS  = ', nconfig
        

     ALLOCATE(cx(ndim,ndim))

     
     cx    = 0.0D+00
     work  = 0.0D+00
     ntest = 0
     nk    = 0

     DO i = 1, ntot

        ntest  =  ntest + ih(i)
        nk     =  nk    + ih(i)
        nkk    =  0

        DO j = 1, ntot

           nkk           = nkk          + ih(j)
           cx(nk, nkk)   = cx(nk,nkk)   + hx(i,j) * dfloat(ih(i)*ih(j))
           work(nk, nkk) = work(nk,nkk) + bx(i,j) * dfloat(ih(i)*ih(j))

           end do
        end do

!!!.............................

        deallocate(ih)
        

        hx = 0.0D+00
        bx = 0.0D+00


        IF(ntest.NE.nconfig) THEN
           WRITE(*,*) '# h2e wrong ntest should equal ntest'
           WRITE(*,*) '# h2e                             ntest = ', ntest 
           WRITE(*,*) '# h2e                           nconfig = ', nconfig 
        ENDIF


        do i = 1, nconfig
           do j = 1, nconfig

              hx(j, i) = cx(j, i)
              bx(j, i) = work(j, i)

           end do
        end do

        WRITE(*,*) '# h2e::       new matrix built.'
        WRITE(*,*) '# h2e::                          row  = ', nk
        WRITE(*,*) '# h2e::                        column = ', nk



        WRITE(*,*)  '# h2e::                      hx matrix:'
        DO n = 1, nconfig, nconfig - 1

           WRITE(*,'(5e12.5)') (hx(n,m),m = 1, nconfig, nconfig-1)

        END DO

        DO n = nconfig - 2, nconfig

           WRITE(*,'(3e12.5)') (hx(n,m),m = n, nconfig)

        END DO


        WRITE(*,*)  '# h2e::                      bx matrix:'

        DO n = nconfig - 2, nconfig

           WRITE(*,'(3e12.5)') ( bx(n,m), m = n, nconfig)

        END DO

!!!..............................

        !         setting up constraints
        ncol = 0
        WRITE(*,*) '# h2e::              initial values:'
        WRITE(*,*) '# h2e:                        nmx = ',  nmx
        WRITE(*,*) '# h2e:                    nconfig = ', nconfig
        WRITE(*,*) '# h2e:                    nconstr = ',  nocnstr
        WRITE(*,*) '# h2e:                       ncol = ',  ncol
       
        WRITE(*,*) '# h2e:: set configuration space.'

        CALL set_space(nocnstr, ncsmx, noll(1)+1)

        WRITE(*,*) '# h2e:: set constraints. '

        CALL constraints(nmx, ncsmx, nconfig, nocnstr, ncol, ncore)

!!!       deallocate(overlap)

        WRITE(*,*) '# h2e::               final values: '
        WRITE(*,*) '# h2e:                        nmx = ',  nmx
        WRITE(*,*) '# h2e:                    nconfig = ', nconfig
        WRITE(*,*) '# h2e:                    nconstr = ',  nocnstr
        WRITE(*,*) '# h2e:                       ncol = ',  ncol


       
       work = 0.0D+00
       DO i = 1, nocnstr

          work(ix(i):ix(i)+ndi(ic(i))-1,iy(i)) = cvec(1:ndi(ic(i)),i)
          work(iy(i),ix(i):ix(i)+ndi(ic(i))-1) = cvec(1:ndi(ic(i)),i)

       END DO

       DEALLOCATE(ix) 
       DEALLOCATE(iy) 
       DEALLOCATE(ic)
 
       !        write(*,*)'ncol=',ncol
       !        write(*,*)'work'
       !        do n=1,nconfig+ncol
       !          do m=1, nconfig+ncol
       !  if(work(m,n).ne.0.0e0) write(*,'(2I5, e12.4)') m,n, work(m,n)
       !          end do
       !        end do
       

!!!..........................


       !                       continuum energy loop

        nconfig = nconfig + ncol

        WRITE(*,*) '# h2e::                     nconfig = ', nconfig
        WRITE(*,*) '# h2e::                        ncol = ', ncol


!!!
!!!    NOTE :: NH2ESINGLE SHOULD BE ALWAYS = 1 FOR SINGLE IONIZATION DATA
!!!            IF TO BE CORRECTED ROUTINE 'H2EFILE' SHOULD MODIFIED 
!!!            ACCORDINGLY



        CALL H2EFILE(NH2ESINGLE, LO) 
        CALL H2EFILE(NH2EDOUBLE, LO) 


        WRITE(NH2ESINGLE) ne
        WRITE(NH2EDOUBLE) ne

!!! TESTING THE ENERGY LOOP 
!!!..............................


        DO i = 1, ne

           eTEST = e0 + de * dfloat( i - 1)

              nd = 0
              DO j = 1, ncsmx

                 IF(idcs(j).EQ.1) THEN

                    IF(ABS(eTEST - eth(j)) < 1.0D-08) THEN
                       WRITE(*,*) '# h2e:: state energy equal to core (threshold) energy'
                       WRITE(*,*) '# h2e::                          en = ', eTEST
                       WRITE(*,*) '# h2e::                      eth(j) = ', eth(j) , j
                       WRITE(*,*) '# h2e:: STOP.'
                       STOP
                    ENDIF
                 END IF

              ENDDO

        ENDDO

!!!................................


        cx = 0.0D+00
        
        CALL ite_space(ndim, ncsmx)

        DO i = 1, ne

           energy = e0 + de * dfloat(i-1)

           WRITE(*,*) '#'
           WRITE(*,*) '# h2e::                    state  i = ', i
           WRITE(*,*) '# h2e::                          en = ', energy

           nd = 0
           DO j = 1, ncsmx
              
              IF(idcs(j).EQ.1) THEN
                 IF(energy > eth(j)) THEN 

                    nd = nd + 1

                    WRITE(*,*) '# h2e::        threshold index    j = ' , j
                    WRITE(*,*) '# h2e::        thrshold energy  eth = ' , eth(j)

                 END IF
              END IF
           END DO

           WRITE(*,*) ' # h2e::         nof channels    nd =' , nd 



!!           OPEN(19,file='cx.dat') 

           !              matrix 'A'
           DO n = 1, nconfig
              DO m = 1, nconfig

                 cx(m,n) = hx(m,n) - energy * bx(m,n) + work(m,n)

!!                 IF(m.le.10.and.n.le.10) THEN
!!                 WRITE(19,*) m, n, cx(m,n) 
!!                 ENDIF

              END DO
           END DO
         
!!           close(19)  

           WRITE(*,*) '# h2e::                   cx(1,1) = ', cx(1,1), adiag(1), ndi(1)
           WRITE(*,*) '# h2e::       cx(nconfig,nconfig) = ', cx(nconfig,nconfig), adiag(1), ndi(1)

         !            write(*,*) adiag(j), ndi(j)
         !            write(*,*) 'cx()=', cx(nmx,nmx)
         !            cx(nmx,nmx)=cx(nmx,nmx)-adiag(j)

           nmx  = 0
           DO j = 1, nd

              nmx         = nmx + ndi(j)  
              cx(nmx,nmx) = 0.0D+00

           END DO
         
                             !!! Old:                 inverting matrix 'A'
                             !  From here on new
         
           WRITE(*,*) '# h2e::                       ndim = ', ndim
           WRITE(*,*) '# h2e::                    nconfig = ', nconfig


         !  factorize the matrix cx ='A'


         
           !          do j=1,nconfig,61
           !             print*,'Matrix=', cx(j,j)
           !          end do
           !
           !
           !  Solve the linear system of Eqns : A * X = B 
           !
           !
           !  A          :      cx(ndim,nconfig) 
           !  LDA        = ndim
           !  order of A = nconfig
           !
           !
           !  SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           !
           !  -- LAPACK driver routine (version 3.0) --
           !      INTEGER            INFO, LDA, LDB, N, NRHS
           !      INTEGER            IPIV( * )
           !      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )




!!!essl
!!$         call dgef(cx,ndim,nconfig, ipvt)
!!$                                !  solve for state vectors
!!$         nmx = 0; er = 0.0d0
!!$
!!$         do j  = 1,nd
!!$
!!$            nmx       = nmx + ndi(j)
!!$            er(nmx,j) = 1.0d0
!!$
!!$            call dges(cx,ndim,nconfig,ipvt,er(:,j),0)
!!$
!!$         end do
!!$
!!$         deallocate(ipvt)
!!!essl




!!!................. SOLVE FOR STATE VECTORS 


!!! LAPACK/NAG
!!!    cx(ndim, ndim) , LDA = ndim  , er(ndim, nchl) , LDB = ndim
!!!    cx ---> cx(nconfig,nconfig)  , er ---> er(nconfig, nd)


           WRITE(*,*) '# h2e:: factorizing hamiltonian.'

           nmx = 0 
           er = 0.0D+00
           DO j  = 1, nd
              nmx   = nmx + ndi(j)
              er(nmx, j) = 1.0D+00
           ENDDO

           nrhs = nd

           WRITE(*,*) '# h2e:: inverting hamiltonian.'



!           CALL dsytrf('u',nconfig,cx,ndim,ipvt,lapack_work,lwork,info) 
!           CALL dsytrs('u',nconfig,nrhs,cx,ndim,ipvt,er,ndim, info) 

           lwork = nconfig * ILAENV(1, 'DSYTRF', 'u', nconfig, -1, -1, -1)

           ALLOCATE(lapack_work(lwork),ipvt(nconfig))

           CALL dsysv('u',nconfig,nrhs,cx,ndim,ipvt,er,ndim, &
                lapack_work, lwork, info)

           DEALLOCATE(ipvt,lapack_work)

           WRITE(*,*) '# h2e:: hamiltonian inversion done.'

!!!outpout of dg2ec.f programm
!!!
!!!
!!!        Energy                   : energy
!!!        Number of open channels  : nd
!!!        Coefficients of Y(E)     : er(n,j)
!!!
!!!        Final result gives  'nd' independent solutions for the wf Y(E)  
!!!       

!           IF(energy.LT.0.0D+00) THEN

              NH2E = NH2ESINGLE
!           ELSE
            
!              NH2E = NH2EDOUBLE
!           ENDIF
         
            
           WRITE(NH2E) energy, nd
           WRITE(NH2E) ((er(n,j), n = 1, nconfig - ncol ), j = 1, nd)


           WRITE(*,*) '# h2e::                         en = ', energy
           WRITE(*,*) '# h2e::                         nd = ', nd    

          
           DO J = 1, ND 
              WRITE(*,*) '# h2e::                  channel j = ', j    
              WRITE(*,*) (er(n,j), n = 1, 10)
              WRITE(*,*) '..........'
           ENDDO
              WRITE(*,*) '#'

!!!           WRITE(*,*) energy, nd
!!!           WRITE(*,*) ((er(n,j), n = 1, nconfig - ncol ), j = 1, nd)
!!           WRITE(*,'(10e11.3)') (er(n,1), n = nmx + 1,nmx + ndi(j))

!!!write(*,'(10e11.3)') (er(n,1), n = 1, nconfig - ncol, nd)

        END DO


!!!.........................
        CLOSE(NH2ESINGLE)
        CLOSE(NH2EDOUBLE)
      
        DEALLOCATE(er)
        DEALLOCATE(vec)
        DEALLOCATE(adiag)
        DEALLOCATE(work)
        DEALLOCATE(hx)
        DEALLOCATE(bx)
        DEALLOCATE(cx)
        DEALLOCATE(eth)
        

!!!........................

    END PROGRAM H2E

!!!###################################################################
!!!   setting up constraints
    SUBROUTINE constraints( nbs, ncsmx, nconfig, ndim, ncol, ncore)

      USE cnstr_data
      USE ioroutines
      IMPLICIT NONE
      
      INTEGER nbs, nconfig, noc, ncol, n0, k, lang, nmx, i, j, l, m, n
      INTEGER ncsmx, ndim, ncore
      INTEGER NOUT, NBMX, ND1E
      REAL(DPK) dr0, rmax, energy
      
!!!      CHARACTER(len=16), ALLOCATABLE, DIMENSION(:) :: fc, hfile
     
      INTEGER, ALLOCATABLE, DIMENSION(:)     :: nc, lc, jc, ncs
      REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: coe
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: bmat

!!!..........................

      NOUT = 16
      ND1E = 1
      NBMX = 8      
      
!!!..........................

      WRITE(*,*) '# h2e::constraints:              ncore = ',  ncore
      WRITE(*,*) '# h2e::constraints:              ndim  = ',  ndim
      WRITE(*,*) '# h2e::constraints:              ncsmx = ',  ncsmx

      allocate( ncs(ncore) )
      allocate( nc(ndim)   )
      allocate( lc(ndim)   )
      allocate( jc(ndim)   )


!..............................

      DO i = 1, ncore
         ncs(i) = 0
      END DO

!!!............................

      ! set here new number of constraints

      noc = 0
      DO i = 1, ncsmx
         
         ncs( lhf(i) + 1 ) = nhf(i)

!!!         WRITE(*,*) 'i,idcs(i) = ',  i, idcs(i)

         SELECT CASE( idcs(i) )

         CASE(1)

WRITE(*,*) '# h2e::constraints:  noc, ncs( ll(i) + 1 )  = ', noc, ncs( ll(i) + 1 )

         DO j = noc + 1, noc + ncs( ll(i) + 1 )

            nc(j) = j - noc
            lc(j) = ll(i)
            ic(j) = i
            jc(j) = j

!         WRITE(*,*) 'j, lc(j), ll(i) ',j, lc(j), ll(i)

      END DO

         !              ncs(lhf(i)+1) = nhf(i) 

         noc = noc + ncs( ll(i) + 1 )

      END SELECT

   END DO
      

   WRITE(*, *) '# h2e::constraints:               noc = ', noc
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        write(*,*) (fc(i),i=1,noc)
!!!        write(*,*) (nc(i),i=1,noc)
!!!        write(*,*) (lc(i),i=1,noc)
!!!        write(*,*) (ic(i),i=1,noc)
!!!        write(*,*) (jc(i),i=1,noc)
!!!        stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    noc  = number of constraints
!!!    fc() = file names of input coef., (xxcoe.??)
!!!    nc() = main quantum numbers for constr.-orbitals
!!!    lc() = ang. quantum numbers for constr.-orbitals
!!!    ic() = the index of the channel constrained
!!!    jc() = the collum index of the constraints
!!!

      ALLOCATE( bmat(nbs, nbs) )
      ALLOCATE( coe(nbs) )


      ncol = 0
      ix   = 0
      iy   = 0
      cvec = 0.0D+00
      bmat = 0.0D+00
      coe  = 0.0D+00

! read one-electron coefficients

!!!.................................


      DO l = 1, noc

         WRITE(*,*) '# h2e::constraints:         l, lc(l) = ', l, lc(l)

         IF(jc(l).GT.ncol) ncol = jc(l)


         CALL D1EFILE(ND1E, lc(l)) 

         READ(1) dr0, rmax, n0, k, lang

         IF(lang.NE.lc(l)) WRITE(*,*) '# h2e::constraints: incosistency in angular momentum.'

         DO i = 1, n0 - 1 - nc(l)

            READ(ND1E) energy, (coe(j),j = 1, n0 - 2)

         END DO

         CLOSE(ND1E)
          
!!!..........................


         !!! read overlap matrices

         CALL BMXFILE(NBMX, lang) 
         
         READ(NBMX) nmx
                       !          write(*,*) 'nmx=',nmx, 'l='
         IF(nmx.NE.nbs) WRITE(*, *) '# h2e::constraints:  incosistency in nmx.'

         !            write(*,'(10e8.2)') (bmat(i,j), j=i,nmx)

         DO i = 1, nmx

            READ(NBMX) (bmat(i,j),j=i,nmx)

            bmat(i:nmx,i)=bmat(i,i:nmx)

         END DO
         
         CLOSE(NBMX)       
          
!!!............................

         DO j = 1, ndi(ic(l))
            DO i = 1, nmx - 1                
               cvec(j,l) = cvec(j,l) + bmat(j,i) * coe(i)
            END DO
         END DO

!!!.............................

         DO m = 1, ic(l) - 1
 
            ix(l) = ix(l) + ndi(m)

         END DO

         ix(l) = ix(l) + 1
          iy(l) = nconfig + jc(l)

!!!          write(*,*) 'cvec', l, ix(l), iy(l)
!!!          write(*,'(10e8.2)') (cvec(j,l),j=1,nmx)

       END DO
       
       DEALLOCATE(nc)
       DEALLOCATE(lc)
       DEALLOCATE(jc)
       DEALLOCATE(coe)
       DEALLOCATE(bmat)
       
       ndim = noc
       
     END SUBROUTINE constraints

!!!####################################################

     SUBROUTINE iterate( nx, nd)   
       USE ite_data
       IMPLICIT NONE
       INTEGER nx, nd, i, j, n, m
       REAL(DPK) fnorm
       
       fnorm = 1.0D+000
       DO i = 1, nd
          
          er(:,i) = vec(:,i)

          DO j = 1, 1
             
             finit(:nx) = er(:nx,i)

                  !              er(:nx,i)= matmul(cx(:nx,:nx),finit(:nx))

             DO n = 1,nx
                er(n,i) = 0.0D+000

                DO m = 1, nx

                   er(n,i) = er(n,i) + cx(m,n) * finit(m)

                END DO
            end do
            !            finit(:nx)=er(:nx,i)


            fnorm = 0.0D+00
 
            DO n = 1, nx
               DO m =1, nx
                  
                  fnorm = fnorm + er(m,i)  *bx(m,n) * er(n,i)

               END DO
            END DO

            WRITE(*,*) '# h2e::iterate:         iteration j =', j, dsqrt(fnorm)

            er(:nx,i) = er(:nx,i) / dsqrt(fnorm)

            WRITE(*,'(10e11.3)') (er(n,i), n=1, nx, 42)

         END DO
      END DO

    END SUBROUTINE iterate
!!!###################################################################
!!!EOF



