!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE param

  USE PRECISION, ONLY : DPK
  IMPLICIT NONE
  PUBLIC


  INTEGER, PARAMETER :: kx = 15, nl = 15, np = 2000
  INTEGER ntc, nsp, id


!..........  h1e.inp file   ...........

  REAL(dpk) znuc, rmin, rmax, xi
  INTEGER   no, idbsp, nrp 
  REAL(dpk) redmass 
  INTEGER   nos,  nop,     nod
  INTEGER   isl,  nw 
  INTEGER   lmin, lcore,   lmax, itrmx
  REAL(dpk), DIMENSION(:),ALLOCATABLE  :: alp1, r01, zk 
  INTEGER NB, KB, IB
  CHARACTER(LEN=3) METHOD 
  INTEGER SELECT
  CHARACTER(LEN=3) GAUGE
  CHARACTER(LEN=2) SPECTRUM_TYPE
  CHARACTER(LEN=3) DISCRETE
  REAL(dpk), DIMENSION(nl)  :: rs
!  REAL(dpk) crt, di, df
  REAL(dpk) en_1e, de_1e
  INTEGER   n_continuum_states, n_bound_states
  REAL(dpk), DIMENSION(:,:), ALLOCATABLE :: rpw



CONTAINS

!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  SUBROUTINE input

    IMPLICIT NONE
    INTEGER i

!!!...............



    OPEN(14,file ='inp/h1e.inp',status='old')

    READ(14, *) nb, kb              ! number of B-splines, order of B-splines 
    READ(14, *) rmin, rmax, no      ! 
    READ(14, *) idbsp               ! selects the grid 0 = sin, 1=exp
    READ(14, *) xi, nrp             ! exp-like
    READ(14, *) znuc, redmass       ! Atomic number, reduced mass
    READ(14, *) nos, nop, nod       ! number of core-orbitals 
    READ(14, *) isl, itrmx          ! HF parameter
    READ(14, *) lmin, lmax, lcore   ! lmin, lcore, lmax

    nw = lmax - lmin + 1            ! number of symmetries

!    ALLOCATE(  RS(NW) )
    ALLOCATE( ALP1(NW) )
    ALLOCATE( R01(NW)  )
    ALLOCATE( ZK(NW)   )

    READ(14, *) ( rs(i),  i = 1, nw)   ! sine-like first non-zero knot points
    READ(14, *) ( alp1(i),i = 1, nw)   ! core-polarizability 
    READ(14, *) ( r01(i), i = 1, nw)   ! cut-off radius for core-polarization or model pot.
    READ(14, *) ( zk(i),  i = 1, nos + nop + nod) ! effective radius for the core shells

    ! inverse iteration(INV), linear solution(LIN)
    
    READ(14, *) METHOD, DISCRETE, select 

    ! initial c. e., e. step, # c. states, b. states
    
    READ(14, *) en_1e, de_1e, n_continuum_states, n_bound_states 
    READ(14, *) GAUGE, SPECTRUM_TYPE

    CLOSE(14)

    ALLOCATE( RPW(NW, NO) )

    nsp = nos + nop
    ntc = nsp + nod

    WRITE(*, *)  '#          1-E FREE BOUNDARY CODE          '
    WRITE(*, *)  '#'
    WRITE(*, *)  '#        NOF  BSPLINES                NB = ', NB
    WRITE(*, *)  '#             ORDER BSPLINES          KB = ', KB
    WRITE(*, *)  '#             BOX   RADIUS             R = ', RMAX
    WRITE(*, *)  '#             ATOMIC NUMBER            Z = ', ZNUC
    WRITE(*, *)  '#        NOF  PARTIAL WAVES           NL = ', NW
    WRITE(*, *)  '#             INITIAL ENERGY         E_0 = ', EN_1E
    WRITE(*, *)  '#             ENERGY  STEP            DE = ', DE_1E
    WRITE(*, *)  '#             LINEAR METHOD  B+C  STATES = ', SELECT
    WRITE(*, *)  '#        NOF  BOUND  STATES          NBS = ', N_BOUND_STATES
    WRITE(*, *)  '#        NOF  CONTINUUM STATES       NCS =',  N_CONTINUUM_STATES
    WRITE(*, *)  '#             DIPOLE M.ELEMENTS    GAUGE = ', GAUGE
    WRITE(*, *)  '#                          SPECTRUM TYPE = ', SPECTRUM_TYPE
    WRITE(*, *)  '#                                    ISL = ', ISL
    WRITE(*, *)  '#                              PRECISION = ', DPK

    IF((ISL.NE.1).AND.(ITRMX.NE.1)) THEN

       WRITE(*,*)'#             HARTREE-FOCK  CALCULATION   '
       WRITE(*,*)'#'
       WRITE(*,*)'#                                  LCORE = ', LCORE
       WRITE(*,*)'#        NOF  S,P,D ORBITALS             = '
       WRITE(*,*)'#                                    NOS = ', NOS
       WRITE(*,*)'#                                    NOP = ', NOP
       WRITE(*,*)'#                                    NOD = ', NOD
       WRITE(*,*)'#        NOF  ITERATIONS           ITER  = ', ITRMX

    ELSE IF((ISL.NE.1).AND.(ITRMX.EQ.1)) THEN

       WRITE(*,*)'#        EXISTING CORE COEFF   USED         '
    ELSE 
       
       IF(REDMASS.EQ.1) THEN
       
          WRITE(*,*)'#  HYDROGENIC RADIAL  FUNCTIONS   CALCULATED'
       ELSE IF(REDMASS.EQ.2) THEN

          WRITE(*,*)'#  POSITRONIUMS  RADIAL FUNCTIONS CALCULATED'
       ELSE
 
          WRITE(*,*)'#    NOT ACCEPTABLE VALUE FOR PARAM REDMASS '
          WRITE(*,*)'#    EXITING...                             '
          STOP

       ENDIF

    ENDIF

    IF((ALP1(LMIN).EQ.0.D+00)) THEN

       WRITE(*,*)'#          NO CORE -POLARIZATION '
    ELSE
       WRITE(*,*)'#           CORE - POLARIZATION IS INCLUDED'
    ENDIF

       WRITE(*,*)'################################################################'

  END SUBROUTINE input

END MODULE param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



