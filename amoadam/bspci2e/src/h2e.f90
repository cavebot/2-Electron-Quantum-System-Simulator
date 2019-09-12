!
!     h2e.f90:: single-processor version
!
! mod_h2e.f90::     cnstr_data:
!                          sub:  set_space: set_overlap_file(1), set_config(2), set_cnstr_array(3)
!         (not used)  ite_data:
!                          sub: set_ite_data   

  PROGRAM H2E
    !
    USE PRECISION, only: dpk
    USE utils,     only: datafile
    USE cnstr_data
    USE ite_data 
    USE ioroutines
    USE bs_frb_2e, only: read_target_energies
    !
    IMPLICIT NONE
    !
    INTEGER::  j, i, n, m, nmx, nmx2,icf, nx, ncol, nd, ne, nk, nkk,nrhs
    INTEGER::  ntot, ncsmx, lo, ls, nconfig, ntest, nchl, inpl
    INTEGER::  ncore, ndmax
    INTEGER::  nocnstr, ndim, ns, nl, nbsp
    INTEGER::  info
    INTEGER::  ifile, nout
    INTEGER::  NCFG, N1E, NR12, NHMXFILE, NBMXFILE
    INTEGER::  NH2E
    REAL(DPK)::  energy, e_1st, de, e_Last, etest
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ih, ipvt
    INTEGER :: lwork, ilaenv
    REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: lapack_work
    REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: adiag, eth
    REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: work, hx!, ehf
    REAL(DPK), DIMENSION(:,:),  POINTER    :: ehf
    CHARACTER(LEN = 100)                   :: ARGV
    !exe!


    ifile      = 8
    NOUT       = 16

    
    NH2E       = 1                     !!! DO NOT CHANGE
    NR12       = 3 
    NHMXFILE   = 7
    NBMXFILE   = 8
    NCFG       = 15
    N1E        = 18  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*, '(a60,i10)') ' h2e start.'

      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL         
      CALL GETARG(2, ARGV)
      READ(ARGV,*) NE         


    OPEN(nout, FILE='out/h2e.out')


    OPEN(ifile, file='inp/h2e.inp')
    READ(ifile,*) e_1st                 ! in Rydberg
    READ(ifile,*) de                    ! in Rydberg
    READ(ifile,*) nocnstr               ! nof constraints
    READ(ifile,*) ncore                 !
    CLOSE(ifile)

    ! read configuration file


    lo = inpl
    CALL getnumberofchannels(ncfg, lo, ndmax, nbsp)
    CALL set_space(NDMAX, NCHL)

    WRITE(*, '(a60,i10)') '                           ndmax = ', NDMAX

    CALL cxfin(ncfg, nhf,lhf,ll,nmin,nmax,noll,is,lo,ls,ncsmx,ntot,ndi,idcs)


    ns = MAXVAL(nmax)
    nl = max0(MAXVAL(lhf), MAXVAL(ll)) + 1

    
    CALL read_target_energies(ehf)                 ! read 1e energies (in Ryd)
        
    WRITE(*,'(a60,i10)') '           nof target    waves = ', SIZE(ehf,dim=1)
    WRITE(*,'(a60,i10)') '           nof target energies = ', SIZE(ehf,dim=2)
    !read 1-e target energies



!    e_1st = 2.0D+00 * e_1st            !go on with Ryd
!    de = 2.0D+00 * de


    !ndi(i)  number of cfg in channel i
    !ndim total number of cfg

    nconfig = 0
    DO i = 1, ncsmx

       n = is(i) - 1

       DO j = 1, ndi(i)   
          n = n + 1
          nconfig = nconfig + 1
       END DO
    END DO
    
    !
    !

    ndim = MAX(ntot, nconfig + nocnstr)


    WRITE(*, '(a60,i10)' ) 'ntot  = ', ntot
    WRITE(*, '(a60,i10)' ) 'nconfig  = ', nconfig
    WRITE(*, '(a60,i10)' ) 'nconstr  = ', nocnstr
    WRITE(*, '(a60,i10)' ) '(nconfig+noconstr)  = ndim  = ', ndim

        

        !allocations for h(2e) / B(2e)

    ALLOCATE( adiag(ncsmx) )
    ALLOCATE( hx(ndim, ndim) )
    ALLOCATE( bx(ndim, ndim) )
    ALLOCATE( eth(ncsmx) )
        
    WRITE(*,*) 'A, H, B allocation done. '
    
    hx    = 0.0D+00
    bx    = 0.0D+00
    adiag = 0.0D+00
    
    DO icf = 1, ncsmx
       
       SELECT CASE(idcs(icf))
              
       CASE(1) ! free boundaries
              
          n = ll(icf) + 1
          
          !read  h(1e) /B(1e)
          
          CALL datafile(nhmxfile, ll(icf), 'hb-')
          CALL datafile(nbmxfile, ll(icf), 'bb-')

              !              CALL HMXFILE( NHMXFILE, ll(icf) ) 
              !              CALL BMXFILE( NBMXFILE, ll(icf) ) 
              
              
          READ(NHMXFILE) nmx, nmx2
          READ(NBMXFILE) nmx, nmx2
              
          DO i = 1, noll(icf)
              
             n   = is(icf) + i - 1
             nx  = is(icf) + noll(icf) - 1
             
             READ(NHMXFILE) ( hx(n,j), j = n, nx)
             READ(NBMXFILE) ( bx(n,j), j = n, nx)
             
                 !                           build H for uncoupled channels 
                 !2.0_dpk because it reads in a.u.
                 !
             hx(n, n:nx) = 2.0_dpk * hx(n,n:nx) + bx(n,n:nx) * ehf( lhf(icf) + 1, nhf(icf) )
             bx(n:nx, n) = bx(n,n:nx)
             
          END DO
          
          adiag(icf) = hx(nx, nx)
              
       CASE(0)  ! fxd boundaries
          
          DO i = 1, noll(icf)
             
             n       = is(icf) + i - 1
             nk      = i + nmin(icf) - 1
             hx(n,n) = ehf(lhf(icf) + 1, nhf(icf)) + ehf(ll(icf)+1, nk)
             bx(n,n) = 1.0D+00
              
          END DO
          
       END SELECT
       
       CLOSE(NHMXFILE)
       CLOSE(NBMXFILE)
       
       eth(icf) = ehf( lhf(icf) + 1, nhf(icf))
    END DO



!      add V_12 to total Hamiltonian matrix 
!
!       H = H(1) + H(2) + H_12 


!!    ndim < dimension of hxh matrix (calculated from r12b.f) program
!!    or total number of configuration  ncf in  hxh.din should be more 
!!    than total number of configurations in dg2ec.din + no of constraints.
!!
!!       ncf(hxh.din) >= ncf(dg2ec.din) + ncostrnt




     ALLOCATE(work(ndim,ndim)) 

     work = 0.0D+00


! read V_12

     CALL hfile(nr12,"dat","r12b","bin",lo)


     DO i = 1, ntot
 
        READ(NR12) ( work(i,j),j = i, ntot)

        hx(i, :)         = hx(i, : )      +  work(i, : )
        hx(i:ntot, i )   = hx(i, i:ntot )
     END DO
     CLOSE(NR12)

     
     ALLOCATE(ih(ndim))

     ih  = 0
     DO i = 1, ncsmx
        
        n = is(i) - 1
           
        DO j = 1, ndi(i)   

           n     =  n + 1
           ih(n) =  1
        END DO
     END DO

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
           cx(nk, nkk)   = cx(nk,nkk)   + hx(i,j) * dble(ih(i)*ih(j))
           work(nk, nkk) = work(nk,nkk) + bx(i,j) * dble(ih(i)*ih(j))

           end do
        end do

        DEALLOCATE(ih)


        !
       
        hx = 0.0_dpk
        bx = 0.0_dpk
        DO i = 1, nconfig
           do j = 1, nconfig
              hx(j, i) = cx(j, i)
              bx(j, i) = work(j, i)
           end do
        end do


        ncol = 0        
        CALL set_space(nocnstr, ncsmx, noll(1)+1)
        CALL constraints(nmx, ncsmx, nconfig, nocnstr, ncol, ncore)

!!!       deallocate(overlap)
       
       work = 0.0D+00
       DO i = 1, nocnstr

          work(ix(i):ix(i)+ndi(ic(i))-1,iy(i)) = cvec(1:ndi(ic(i)),i)
          work(iy(i),ix(i):ix(i)+ndi(ic(i))-1) = cvec(1:ndi(ic(i)),i)

       END DO

       DEALLOCATE(ix) 
       DEALLOCATE(iy) 
       DEALLOCATE(ic)
       !                       continuum energy loop

        nconfig = nconfig + ncol
        test_en_loop:DO i = 1, ne
           eTEST = e_1st + de * DBLE( i - 1)
              nd = 0
              DO j = 1, ncsmx
                 IF(idcs(j).EQ.1) THEN
                    IF(ABS(eTEST - eth(j)) < 1.0D-08) THEN
                       WRITE(*,*) 'state energy equal to core (threshold) energy'
                       WRITE(*,*) '     en = ', eTEST
                       WRITE(*,*) ' eth(j) = ', eth(j) , j
                       WRITE(*,*) ' stop'
                       STOP
                    ENDIF
                 END IF
              ENDDO
           ENDDO test_en_loop


!           nbound = 0
!           get_bound_states: DO i = 1,  ne
!              energy = ehf(2,i)
!              IF(energy.LT.ehf(1,8)) THEN 
!                 nbound = nbound + 1
!                 CYCLE
!              ENDIF
!           ENDDO 
              
           !

           CALL hfile(nh2e,"dat","h2e","bin",lo)

           WRITE(NH2E) ne

           cx = 0.0_dpk
        
           CALL ite_space(ndim, ncsmx)


           e_Last = e_1st + ne * de 

           WRITE(*,'(a60)')'&h2e:: solve for the chosen energy grid:'
           WRITE(*, '(a60,E15.5)') 'e_1st  = ', e_1st
           WRITE(*, '(a60,E15.5)') 'de  = ', de
           WRITE(*, '(a60,i10  )') 'ne  = ', ne
           WRITE(*, '(a60,E15.5)') 'e_Last  = ', e_Last

           find_state_at_energy: DO i = 1, ne

              energy = e_1st + de * DBLE(i-1)               !energy grid

!           find_state_at_energy: DO i = 1,  ne !ns!SIZE(ehf(lhf(1)=1,nhf(1))
!              energy = ehf(2,7+i)
!              IF(energy.LT.ehf(1,8)) THEN                 
!                 WRITE(*,*) '     i,en = ', i,energy, ehf(1,8)
!                 CYCLE
!              ENDIF

              WRITE(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
              WRITE(*,'(a60,i10)') 'i = ', i
              WRITE(*,'(a60,E20.10)') ' E = ', energy*0.5_dpk

              nd = 0
              evaluate_nof_open_channels:DO j = 1, ncsmx                 
                 IF(idcs(j).EQ.1) THEN
                    IF(energy > eth(j)) THEN 
                       nd = nd + 1
                    END IF
                 END IF
              END DO evaluate_nof_open_channels

              WRITE(*,'(a60,i10)') 'nof channels    nd(E) =' , nd 

           !              matrix 'A'
              construct_matrix_A:DO n = 1, nconfig
                 DO m = 1, nconfig
                    cx(m,n) = hx(m,n) - energy * bx(m,n) + work(m,n)
                 END DO
              END DO construct_matrix_A
         
              ! set zero the first element in block matrices corresponding to channels
              
              er = 0.0D+00
              nmx  = 0
              set_bc_channels: DO j = 1, nd
                 nmx         = nmx + ndi(j)  
                 cx(nmx,nmx) = 0.0D+00
                 er(nmx, j) = 1.0D+00
              END DO set_bc_channels

!              nmx = 0 
!              er = 0.0D+00
!              DO j  = 1, nd
!                 nmx   = nmx + ndi(j)
!                 er(nmx, j) = 1.0D+00
!              ENDDO

         

              !  solve        (H-E)*C_i = C_0_i


              ! LAPACK/NAG
              !    cx(ndim, ndim) , LDA = ndim  , er(ndim, nchl) , LDB = ndim
              !    cx ---> cx(nconfig,nconfig)  , er ---> er(nconfig, nd)


           nrhs = nd

           WRITE(*,*) 'factorizing/inverting the matrix SEqn .'


           lwork = nconfig * ILAENV(1, 'DSYTRF', 'u', nconfig, -1, -1, -1)

           ALLOCATE(lapack_work(lwork),ipvt(nconfig))
           CALL dsysv('u',nconfig,nrhs,cx,ndim,ipvt,er,ndim, lapack_work, lwork, info)

           DEALLOCATE(ipvt,lapack_work)
           
           WRITE(*,*) 'hamiltonian inversion done.'


           !output:

           !        Energy                   : energy
           !        Number of open channels  : nd
           !        Coefficients of Y(E)     : er(n,j)
           
           !        Final result gives  'nd' independent solutions for the wf Y(E)  


            
           WRITE(NH2E) energy, nd
           WRITE(NH2E) ((er(n,j), n = 1, nconfig - ncol ), j = 1, nd)
          
           DO J = 1, ND 
              WRITE(*,'(a60,i10)') 'channel j = ', j    
              WRITE(*,'(5E15.7)') (er(n,j), n = 1, 5)
           ENDDO

           END DO find_state_at_energy
           CLOSE(NH2E)           


           !deallocations


           DEALLOCATE(ehf)      
           DEALLOCATE(er)
           DEALLOCATE(vec)
           DEALLOCATE(adiag)
           DEALLOCATE(work)
           DEALLOCATE(hx)
           DEALLOCATE(bx)
           DEALLOCATE(cx)
           DEALLOCATE(eth)
           

           ! return

         END PROGRAM H2E
         
!!!##################################################################
!!!   setting up constraints
    SUBROUTINE constraints( nbs, ncsmx, nconfig, ndim, ncol, ncore)
      !mod!
      USE cnstr_data
      USE ioroutines
      USE utils, ONLY : datafile
      !
      IMPLICIT NONE
      !arg!
      INTEGER                                :: nbs
      INTEGER                                :: ncsmx
      INTEGER                                :: nconfig
      INTEGER                                :: ndim
      INTEGER                                :: ncol
      INTEGER                                :: ncore
      !loc!
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: bmat
      REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: coe
      INTEGER,   ALLOCATABLE, DIMENSION(:)   :: nc
      INTEGER,   ALLOCATABLE, DIMENSION(:)   :: lc
      INTEGER,   ALLOCATABLE, DIMENSION(:)   :: jc
      INTEGER,   ALLOCATABLE, DIMENSION(:)   :: ncs       !nof constraints
      REAL(DPK)                              :: energy
      INTEGER                                :: n0
      INTEGER                                :: noc
      !io!
      INTEGER                                :: nout,nbin
      !ind!
      INTEGER                                :: nmx,nmx2  
      INTEGER                                :: i,j,l,m
      !exe!
      !..........................
      
      WRITE(*,  '(a60)'   ) 'h2e:: constraints  in.'     
      WRITE(*, '(a60,i10)') '     noc: nof constraints'
      WRITE(*, '(a60,i10)') '      nc: principal q.n for constraint orbitals,'
      WRITE(*, '(a60,i10)') '      lc: angular   q.n for constraint orbitals,'
      WRITE(*, '(a60,i10)') '      ic: index          of constraint channel'
      WRITE(*, '(a60,i10)') '      jc: index (column) of constraint channel'
      WRITE(*,'(a60,i10)')     'nbs = ',  nbs
      WRITE(*,'(a60,i10)')   'ncsmx = ',  ncsmx
      WRITE(*,'(a60,i10)') 'nconfig = ',  nconfig
      WRITE(*,'(a60,i10)')   'ndim  = ',  ndim
      WRITE(*,'(a60,i10)')    'ncol = ',  ncol
      WRITE(*,'(a60,i10)')   'ncore = ',  ncore

      
      !

      NOUT = 16
      nbin = 1
      
      !

      allocate( ncs(ncore)  )
      allocate(  nc(ndim)   )
      allocate(  lc(ndim)   )
      allocate(  jc(ndim)   )

      !

      ! set here new number of constraints


      ncs = 0 
      noc = 0
      set_constraints: do i = 1, ncsmx
         
         ncs( lhf(i) + 1 ) = nhf(i)

         SELECT CASE( idcs(i) )

         CASE(1)

            DO j = noc + 1, noc + ncs( ll(i) + 1 )
            
               nc(j) = j - noc
               lc(j) = ll(i)
               ic(j) = i
               jc(j) = j
                              
               WRITE(*,'(a60,4i5)') 'j, (nc,lc,jc)  = ', j, nc(j), lc(j), jc(j)
            END DO
            
         !              ncs(lhf(i)+1) = nhf(i) 
            
            noc = noc + ncs( ll(i) + 1 )            
         END SELECT

         WRITE(*,'(a60,4i5)') 'noc_ch, noc, l_ch, i_ch  = ', ncs( ll(i) + 1 ), noc, ll(i), i
         
      END DO set_constraints

      WRITE(*, '(a60,i10)') 'noc = ', noc
      

      !    noc  = number of constraints
      !    fc() = file names of input coef., (xxcoe.??)
      !    nc() = main quantum numbers for constr.-orbitals
      !    lc() = ang. quantum numbers for constr.-orbitals
      !    ic() = the index of the channel constrained
      !    jc() = the collum index of the constraints
!

      ALLOCATE( bmat(nbs, nbs) )
      ALLOCATE( coe(nbs) )


      ncol = 0
      ix   = 0
      iy   = 0
      cvec = 0.0_dpk
      bmat = 0.0_dpk
      coe  = 0.0_dpk

      ! read one-electron coefficients

!!!.................................


      DO l = 1, noc
         
         WRITE(*,'(a60,2i10)') ' noc(i), lc(i) = ', l, lc(l)

         IF(jc(l).GT.ncol) ncol = jc(l)

!         CALL d1efile(nbin, lc(l)) 
!         READ(nbin) n0,nmx2
!         READ(nbin) dr0, rmax, n0, k, lang
!         IF(lang.NE.lc(l)) WRITE(*,*) 'incosistency in angular momentum.'
!         read_target_orbitals_coef: DO i = 1, n0 - 1 - nc(l)
!            READ(nbin) energy, (coe(j), j = 1, n0 - 2)
!         END DO read_target_orbitals_coef
!         CLOSE(nbin)

         CALL datafile(nbin, lc(l), 'h1e-')       ! target orbitals (fxd)
         READ(nbin) n0,nmx2
!         read_target_orbitals_coef: DO i = 1, n0 - 1 - nc(l)
         read_target_orbitals_coef: DO i = 1, nc(l)
            READ(nbin) energy
            READ(nbin)(coe(j), j = 1, n0 - 2)
         END DO read_target_orbitals_coef
         CLOSE(nbin)
!         WRITE(*,*) energy,(coe(j), j = 1, n0 - 2)

!         energy = 2.0_dpk * energy
          
!!!..........................


         !!! read overlap matrices   <i|B|j>
!         CALL datafile(nbmxfile, ll(icf), 'bb-')
!         CALL BMXFILE(nbmx, lang)          
         CALL datafile(nbin, lc(l), 'bb-')        

         READ(nbin) nmx, nmx2
         IF(nmx.NE.nbs) WRITE(*, *) '# constraints:  incosistency in nmx.'
         read_overlap_matrix:DO i = 1, nmx
            READ(nbin) (bmat(i,j), j = i, nmx)
            bmat(i:nmx,i)= bmat(i,i:nmx)
         END DO read_overlap_matrix         
         CLOSE(nbin)
          
         !
         ! note: bmat is a (nb,nb) matrix but with (nb-1, nb-1) non-zero elements
         !       coe  is a (nb)    matrix but with (nb-2) non-zero elements (target orbitals)

         DO j = 1, ndi(ic(l))
            DO i = 1, nmx - 1                
!               cvec(j,l) = cvec(j,l) + bmat(j,i) * coe(i)
               cvec(j,l) = cvec(j,l) + bmat(j,i) * coe(i)
            END DO
         END DO

         ! cvec = matmul(bmat,coe)



         !

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

       WRITE(*,'(a60,i10)') 'h2e:: constraints  out.'
     END SUBROUTINE constraints

!!!####################################################

     SUBROUTINE iterate( nx, nd)
       !
       USE ite_data
       !
       IMPLICIT NONE
       INTEGER nx, nd, i, j, n, m
       REAL(DPK) fnorm
       !EXE!

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

            WRITE(*,*) '# iterate:         iteration j =', j, dsqrt(fnorm)

            er(:nx,i) = er(:nx,i) / dsqrt(fnorm)

            WRITE(*,'(10e11.3)') (er(n,i), n=1, nx, 42)

         END DO
      END DO

    END SUBROUTINE iterate
!!!###################################################################
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


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!EOF



