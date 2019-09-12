!###################################################################
! laan jan/feb 2003 iesl day  : bsp1ef/h1e.f90
!
!
!   Author         :   Lambros Nikolopoulos
!                  :
!!    JAN/FEB 2003 :
!                  :
!   I/0 files      :   inp/h1e.inp 
!                  :   out/h1e.out
!                  :   out/setmat.log
!                  :   out/subsolve.lin.out
!                  :   out/subsolve.inv.out
!                  :
!   I/O data files :   DAT/HD1E-L.LINDAT,   L = 0, 1,  ..., NL
!                  :   DAT/HD1E-L.INV.DAT,  L = 0, 1,  ..., NL
!                  :   DAT/HMX-L.DAT,        L = 0, 1, ..., NL
!                  :   DAT/BMX-L.DAT,        L = 0, 1, ..., NL
!                  :   dat/knot.dat 
!   LIBRARIES      :   NAG/FO2SDF   
!                  :   NAG/F02FHF 
!                  :   NAG/F07BDF, NAG/F07BEF
!                  :
!   SUB/MOD        :  modules.f90, sub_vhf.f, modio.f90
! 
!
!
!    h1e.f90 :
!               1. Calculates energies + coefficients for l = 0, 1, 2,..,nl
!                  for 1-e atomic hydrogenic systems with free boundary 
!                  conditions at the end of the interval [0, R].
!                  Two methods have been employed :
!                  a. inverse iteration              (FO2SDF)
!                  b. solution of linear systems
!
!               2. Calculates and stores the 1-e Hamiltonian on Bsplines 
!                  as well as the B-splines overlap.
!
!                  H_ij = < B_i | -d^2/dr^2 + V(r) | B_j >
!                  B_ij = < B_i | B_j >
!
!  2003.02.05
!
!      to be done: store the matrices in banded form directly
!                  instead of calling full and then transform them into
!                  banded matrices
!
!
! nbs : nof bound states
! ncs : nof continuum states
!


PROGRAM H1E

  USE param
  USE set_grid
  USE one_e_matrix
  USE data, ONLY: write_mx, write_v_mx
  USE ioroutines 

  !..................

  IMPLICIT NONE

  INTEGER  i,j
  INTEGER  l, ns                                   ! nof energy states at symmetry 'l' (nbs + ncs)

  ! I/O

  INTEGER NOUT                            ! I/O FILE S


  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: H, B    ! H_ij, B_ij matrices
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: D
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: C0
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: V_HF_DIRECT, V_HF_EXCHANGE  ! 

  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: P
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: T
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en         ! energy vector i = 1,2,..,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce         ! coefficients for state en(i)
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en_c 
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce_c
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en_b
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce_b
  INTEGER,   ALLOCATABLE, DIMENSION(:)   :: nvalence
  CHARACTER*100 ARGV

  INTERFACE
     SUBROUTINE ENERGY_SPECTRUM(L, EN, MODE)
       USE param, ONLY:dpk
       IMPLICIT NONE
       INTEGER L
       REAL(DPK), DIMENSION(:), INTENT(INOUT) :: en
       CHARACTER(LEN=*)  MODE
     END SUBROUTINE ENERGY_SPECTRUM
     
    SUBROUTINE SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)  
      USE param, ONLY : DPK
      IMPLICIT NONE
      INTEGER L
      REAL(DPK), DIMENSION(:,:),    INTENT(IN) :: H,B
      REAL(DPK), DIMENSION(:,:),    INTENT(IN) :: c0
      REAL(DPK), DIMENSION(:),   INTENT(inout) :: en
      REAL(DPK), DIMENSION(:,:), INTENT(inout) :: ce
    END SUBROUTINE SOLVE_BANDED_LIN_NAG77
    
 END INTERFACE
 !........................... get arguments


  CALL GETARG(1, ARGV)
  READ(ARGV,*)   L 
!........................................


  WRITE(*, '(a2,1X,a42,i3)') "#"," h1e.f90: Rh1e ", L
  WRITE(*, '(a2,1X,a42,i3)') "#"," Input file: inp/h1e.inp"
  WRITE(*, '(a2,1X,a42,i3)') "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


  NOUT         = 16


  !.....................................

! call interfaces

  WRITE(*, '(a45,i3)') " #                      angular symmetry l = ", L

  CALL input
!  CALL cal_rpw
  CALL r_grid


  !.....................................

  ALLOCATE( T(NB + KB) )
  ALLOCATE( NVALENCE( 1 ) )
  ALLOCATE( V_HF_DIRECT(  NB-1, NB-1 ))
  ALLOCATE( V_HF_EXCHANGE(NB-1, NB-1 ))


  NVALENCE    = 1

!.....................................

  OPEN(NOUT,FILE='out/h1e.log')

  
  !  calculate one-e Hamiltonian and B-splines overlap matrix and store !
  !  them in  'amat' and 'bmat' and in dataFiles for later use.         !
  


  ALLOCATE( H(NB-1, NB-1)  )         ! hamiltonian matrix
  ALLOCATE( B(NB-1, NB-1)  )         ! b splines overlap matrix
  ALLOCATE( C0(NB-1, 1)    )


  CALL MKGRID( NB, KB, RMAX, RS, T )


  ! calculate H,B

  CALL SETMAT(NB, L, H, B, T, V_HF_DIRECT, V_HF_EXCHANGE)

  C0         = 0.0D+00
  C0(NB-1, 1) = 1.0D+00
     
  !........   H(i,j) , B(i,j)

  CALL WRITE_MX(L, H, "hmx-")
  CALL WRITE_MX(L, B, "bmx-")
     


! Solve shrodinger 1-e equation



  IF(method.EQ.'l') THEN  ! linear equations method    ( P(R) free )

     ns = nbs + ncs

     ALLOCATE( EN(NS) )
     ALLOCATE( CE(NS, NB-1) )

     CALL ENERGY_SPECTRUM(L, EN, 'OUT')                   ! calculate en(i) and store
     CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)

           
  ELSE IF(method == 'd') THEN ! diagonalization method ( P(R) fxd  (==0) )
        
     ALLOCATE( EN(NB-2) )
     ALLOCATE( CE(NB-2, NB-1) )

     CALL SOLVE_BANDED_DIAG_NAG77(L, H, B, EN, CE) 

     CALL ENERGY_SPECTRUM(L, EN, 'IN')                      ! get en(i) and store

     
  ELSE IF(method == 'dl') THEN ! diagonalization (bound states) + linear (continuum states) method

     ALLOCATE( EN_B(NB-2) )
     ALLOCATE( CE_B(NB-2, NB-1) )

     CALL SOLVE_BANDED_DIAG_NAG77(L, H, B, EN_B, CE_B)     ! solve for en<0 (bound spectrum) 

! get nof of bound states by fxd conditions

     nbs = 0
     DO i = 1, SIZE(en_b)
           
        IF(en_b(i) < 0.0D+00) nbs = nbs + 1
     ENDDO

     WRITE(*,*) ' diag nbs = ', nbs

     IF(spectrum=='fxd') THEN 

! note that here we need only eigenvalues of the fxd proble. To be accomodate it

        ns = SIZE(en_b)

        WRITE(*,*) '         ns = ', ns

        ALLOCATE( EN(NS) )
        ALLOCATE( CE(NS, NB-1) )

        en = en_b

        CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)  ! solve for en>0  (continuum spectrum) 

!??        
        ce(1:nbs,:)    = ce_b(1:nbs,:)

        DEALLOCATE( EN_B )    !en_b
        DEALLOCATE( CE_B )    !ce_b

     ELSE IF (spectrum=='mxd') THEN


        ALLOCATE( EN_C(NCS) )
        ALLOCATE( CE_C(NCS, NB-1) )

        CALL ENERGY_SPECTRUM(L, EN_C, 'INOUT')             ! calculate en_c(i) (do not store)
        CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN_C, CE_C)  !solve for en>0  (continuum spectrum) 
        
        ns = nbs + ncs
        
        WRITE(*,*) "# h1e::                       nbs = ", nbs
        WRITE(*,*) "# h1e::                       ncs = ", ncs

        ALLOCATE( EN(NS) )
        ALLOCATE( CE(NS, NB-1) )
        
! energy spectrum (b + c)

        en(     1 : nbs ) = en_b( 1 : nbs)
        en( nbs+1 :  ns ) = en_c( 1 : ncs)

! coefficients    (b + c)

        ce(    1: nbs,  : ) = ce_b( 1 : nbs, : )
        ce(nbs+1:  ns , : ) = ce_c( 1 : ncs, : )

        DEALLOCATE( EN_B )    !en_b
        DEALLOCATE( CE_B )    !ce_b
        DEALLOCATE( EN_C )    !en_c
        DEALLOCATE( CE_C )    !ce_c

     ENDIF

  CALL ENERGY_SPECTRUM(L, EN, 'IN')  ! get en(i) and store

ENDIF

  CALL WRITE_V_MX(L, EN, CE, "h1e-")


     !........................  DEALLOCATE  ARRAYS FOR L SYMMETRY
  
  DEALLOCATE( NVALENCE)
  DEALLOCATE( T )

  DEALLOCATE( H )
  DEALLOCATE( B )
  DEALLOCATE( C0 )

!! FREE
  DEALLOCATE( EN)
  DEALLOCATE( CE)

  DEALLOCATE( V_HF_DIRECT, V_HF_EXCHANGE )
     
  WRITE(*,*) '# h1e:: Hamiltonian matrix:       h(ij)   = < B(i) | H | B(j) > : dat/hmx.dat'
  WRITE(*,*) '# h1e:: Hamiltonian matrix:       b(ij)   = < B(i) |B(j) >      : dat/bmx.dat'
  WRITE(*,*) '# h1e:: end succesfully.'
  WRITE(*, '(a2,1X,a42,i3)') "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  CLOSE(NOUT)


END PROGRAM h1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF


