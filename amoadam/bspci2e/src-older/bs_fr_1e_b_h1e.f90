!!##################################################################
!!!   f = r^L_(<)/r^{L+1}_(>),   r_(<) = min(x,r/2), r_(>) = max(x,r/2)  

PROGRAM H1EZ

  USE param
  USE set_grid
  USE one_e_matrix
  USE DATA, ONLY: write_mx,write_v_mx
  !..................
  !
  IMPLICIT NONE
  INTEGER  i
  INTEGER  l, ns                                              !
  INTEGER NOUT                                                ! I/O FILE S
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: T                 ! B - splines grid
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: B                 ! total B-splines overlap
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: H0                ! total diagonal (in l,l') hamiltonian 
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: V                 ! coupling l's potential
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en, en_b, en_c    ! energy vectors i = 1,2,..,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce, ce_b, ce_c    ! coefficients for state en(i)
!
  CHARACTER*100 ARGV
!  EXTERNAL f

! executable statements !

  NOUT = 16
  OPEN(NOUT,FILE='out/h1e.log')
!
!
  CALL input                        !get parameters (parameter.f90)
!
!
  CALL GETARG(1, ARGV)              !get arguments
  READ(ARGV,*)   l0                 ! l = 0,2,4... or 1,3,5... 
  CALL GETARG(2, ARGV)              !get arguments
  READ(ARGV,*)   nl 
!....
!

  
  WRITE(*, '(a2,1X,a42,2i3,1X,E15.2)') "#"," hz1e.f90: Rh1e l0, nl, b_f : ", l0, NL, B_F
  WRITE(*, '(a2,1X,a42,i3)')    "#"," Input file: inp/hz1e.inp"
  WRITE(*, '(a2,1X,a42,i3)')    "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  WRITE(*, '(a2,1X,a42,i3)')    "#","           partial wave =", l0
  WRITE(*, '(a2,1X,a42,i3)')    "#"," nof partial waves   nl =", nl
  WRITE(*, '(a2,1X,a42,E15.2)') "#","     magnetic field b_f =", b_f

!
!..    allocations!
!
!

  ALLOCATE(          T(NB + KB) )          !   grid knot points 
  ALLOCATE(  B( NL * (NB-1), NL*(NB-1)) )  ! total hamiltonian
  ALLOCATE( H0( NL * (NB-1), NL*(NB-1)) )  ! total hamiltonian
  ALLOCATE(  V( NL * (NB-1), NL*(NB-1)) )  ! total hamiltonian

!
!...
!
!

  CALL R_GRID                               ! select knot-sequence
  CALL MKGRID(  NB, KB, RMAX, RS, T )       ! make grid 

!
!
! Construct matrix representations of operators on B-spline basis  
!
!      H = H_0 + V       (V: coupling part)
!
! H_0 = -1/2 d^2/dr^2 + l*(l+1)/2r^2 + Za/r 
!
! note: thus the 1-e radial Hamiltonian includes the Coulombing potential of the nuclei.
!  


  CALL  MAKE_B( t, b  )                      ! construct b    
  CALL  MAKE_H0(t, h0 )                      ! construct h0

  IF(code=='bs_fx') THEN 
     h0 = 2.0d+00*h0
  ENDIF

  IF(nl.EQ.1) THEN 
     CALL WRITE_MX(L0, H0, "hmx-")
     CALL WRITE_MX(L0,  B, "bmx-")
  ENDIF


! construct potential v(r)

  v = 0.0D+00                            ! hydrogen (default)

  IF( problem == '1eb') THEN
     CALL  MAKE_V_MAGNETIC( t, v  )      ! Hydrogen in static magnetic field
!     WRITE(*,*) "stop",v
!     stop
  ELSE IF( problem == 'h2p') THEN
     CALL  MAKE_V_MOLECULAR( t, v)       !molecular hydrogen ion (h2+)
  ELSE IF( problem == '2e') THEN

     CALL  MAKE_V_2e( t, v)              !helium
  ENDIF

!
! solve Schrodinger Equation!  
!
!                      H * C = E*B*C  
!      or              A(E)* C = C_0,        A(E) = H - E * B 
!


  IF(method.EQ.'l') THEN              ! linear equations method    ( P(R) free )

     ns = nbs + ncs

     ALLOCATE( EN(NS) )
     ALLOCATE( CE(NS, NL*(NB-1)) )

     CALL  ENERGY_SPECTRUM( L0, EN, 'OUT')             ! calculate en(i) 
     CALL LINEARIZE_F07BEF( H0,  B, EN, CE)            ! solve for en>0
   
  ELSE IF(method == 'd') THEN       ! diagonalization method ( P(R) fxd  (==0) )
        
     ALLOCATE( EN( NL*(NB-2) ))
     ALLOCATE( CE( NL*(NB-2), NL*(NB-1)) )

!     CALL DIAGONALIZE_F02FHF( H, B )                 !get en(i) and C(i)

     CALL DIAGONALIZE_DSBGV( B, H0, V, EN, CE)        !get en(i) and C(i)
     CALL   ENERGY_SPECTRUM( L0, EN, 'IN')             


  ELSE IF(method == 'dl') THEN ! diagonalization (bound states) + linear (continuum states) method


     ALLOCATE( EN_B(NL*(NB-2) ))
     ALLOCATE( CE_B(NB-2, NL*(NB-1)) )

     CALL DIAGONALIZE_DSBGV( B, H0, V, EN_B, CE_B)          !get E(n) and C(n)

     nbs = 0
     DO i = 1, SIZE(en_b)                            ! get nof of b. states by fxd conditions
        IF(en_b(i) < 0.0D+00) nbs = nbs + 1
     ENDDO

     WRITE(*,*) '# solve_banded_diag       nbs = ', nbs

     IF(spectrum=='fxd') THEN  ! note that here we need only eigenvalues of the fxd problem. 
                               ! calculate coe for E_n > 0 by linear-equation solver 
        ns = SIZE(en_b)        ! (results should be the same as in 'd' case

        WRITE(*,*) '         ns = ', ns

        ALLOCATE( EN(NS) )
        ALLOCATE( CE(NS, NL*(NB-1)) )

        en = en_b

        CALL LINEARIZE_F07BEF(H0, B, EN, CE)       ! solve for E(n) > 0  (P(R) = 0) 
                                                   ! and get   C(n) by LE solver 
        ce(1:nbs,:)    = ce_b(1:nbs,:)

        DEALLOCATE( EN_B )    !en_b
        DEALLOCATE( CE_B )    !ce_b

     ELSE IF (spectrum=='mxd') THEN         

        ALLOCATE( EN_C(NCS) )
        ALLOCATE( CE_C(NCS, NL*(NB-1)) )

        CALL ENERGY_SPECTRUM(L0, EN_C, 'INOUT')         !calculate en_c(i) (do not store)
        CALL LINEARIZE_F07BEF(H0, B, EN_C, CE_C)        !solve for E > 0
                                                        !and get  C(E) by LE solver
        ns = nbs + ncs
        
        WRITE(*,*) "# h1e::                       nbs = ", nbs
        WRITE(*,*) "# h1e::                       ncs = ", ncs

        ALLOCATE( EN(NS) )
        ALLOCATE( CE(NS, NL*(NB-1)) )
       
!... assemble  bound and continuum results for E,C(E)

        en(     1 : nbs ) = en_b( 1 : nbs)           ! energy spectrum (b + c)
        en( nbs+1 :  ns ) = en_c( 1 : ncs)
        !
        ce(    1: nbs,  : ) = ce_b( 1 : nbs, : )     ! coefficients    (b + c)
        ce(nbs+1:  ns , : ) = ce_c( 1 : ncs, : )

!...
        DEALLOCATE( EN_B )    !en_b
        DEALLOCATE( CE_B )    !ce_b
        DEALLOCATE( EN_C )    !en_c
        DEALLOCATE( CE_C )    !ce_c

     ENDIF

     CALL ENERGY_SPECTRUM(L0, EN, 'IN')         ! get en(i) and store

  ENDIF

!  DATA_FILE = "h1e-"
!  IF(code == 'bs_fx')   DATA_FILE = "d1e-"      ! P(R) == 0 
!  CALL WRITE_V_MX(L0, EN, CE, DATA_FILE)


  IF(code == 'bs_fx') THEN                 ! P(R) == 0 
     CALL WRITE_V_MX(L0, EN, CE, "d1e-")
  ELSE
     CALL WRITE_V_MX(L0, EN, CE, "h1e-")
  ENDIF


  !........................  DEALLOCATE  ARRAYS FOR L SYMMETRY
  
  DEALLOCATE( T )
  DEALLOCATE( H0,B,V)
  DEALLOCATE( EN)
  DEALLOCATE( CE)


  WRITE(*,*) '# h1e:: Hamiltonian matrix:       h(ij)   = < B(i) | H | B(j) > : dat/hmx.dat'
  WRITE(*,*) '# h1e:: Hamiltonian matrix:       b(ij)   = < B(i) |B(j) >      : dat/bmx.dat'
!  WRITE(*,*) '# h1e:: energies+coefficients in: ',DATA_FILE

  WRITE(*,*) '# h1e:: end succesfully.'
  WRITE(*, '(a2,1X,a42,i3)') "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  CLOSE(NOUT)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (he atomic problem) 
SUBROUTINE MAKE_V_2E(T, V)

REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v

WRITE(*,*) ' h1e:: need to be implemented 31.03.2006'

END SUBROUTINE MAKE_V_2E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (h2+ molecular problem)
SUBROUTINE MAKE_V_MOLECULAR(T, V)

  USE ang,        only : cleb
  USE utils,      only : print_mx
  USE units,      only : m_pi
  USE one_e_matrix, only:f
!................................
!

  IMPLICIT NONE

!
!declarations!
!

  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v
! locals !
  REAL(DPK), DIMENSION(nb-1,nb-1,nl) :: mx_r_lm
  REAL(DPK), DIMENSION(nl,nl,nl)     :: f_gaunt
  INTEGER i,j,il,jl
  INTEGER dim_h_1, dim_h_2
  INTEGER la, lb, lm 
  REAL(DPK) lab
  REAL(DPK) vlm
!  EXTERNAL f

!executable part!....
  WRITE(*,*) '# h1e::make_v_molecular:       call.'
!
!
  dim_h_1 = nl* SIZE(mx_r_lm, 1)
  dim_h_2 = nl* SIZE(mx_r_lm, 1)
!
!...
!

  ! make gaunt factor 

  f_gaunt = 0.0D+00

  la = l0
  DO il = 1, nl

     lb = l0

     DO jl = 1, nl

        
        lab  = (2.0D+00 * DBLE(la) + 1) * (2.0d+00 * DBLE(lb) + 1)
!
        lm = 0
        DO l = 1, nl     ! sum over L
           
           f_gaunt(il,jl,l) = SQRT(lab)                                      &
                            * cleb( 2*la,     0, 2*lm, 0, 2*lb,   0   )      &
                            * cleb( 2*la, 2*m_l, 2*lm, 0, 2*lb, 2*m_l )
           lm = lm + 2
          ENDDO
!
        lb = lb + 2
     ENDDO
!
     la = la + 2
  ENDDO

  !
  ! make molecular potential
  ! mix radial multipoles and angular factors
  !

  lm = 0
  DO l = 1, nl
     
     CALL  mat_bsp(nb, t, mx_r_lm(:,:,l), f, lm, 0 )  !multipole expansion on B-splines
     
     mx_r_lm(:,:,l) = -2*SQRT(4*M_PI/(2*l+1)) * mx_r_lm(:,:,l)

     lm = lm + 2
  ENDDO


  V = 0.0d+00
  
  bsp_i: DO i = 1, nb-1
     bsp_j:   DO j = 1 , nb-1
        DO il = 1, nl 
           DO jl = 1, nl


!              IF( il.NE.jl ) THEN
!                    vlm = 0.0d+00
!                    DO l = 1, nl
!                       vlm = vlm + f_gaunt(il,jl,l) * mx_r_lm(i,j,l) 
!                    ENDDO
!                    V( (i-1)*nl+il, (j-1)*nl+jl) = vlm
!
! check the diagonal and non-diagonal parts L=0 --> 2/r
!

              IF(il.NE.jl) THEN
                 V( (i-1)*nl+il, (j-1)*nl+jl ) =  DOT_PRODUCT(mx_r_lm(i,j,:), f_gaunt(il,jl,:))
              ENDIF

!              ELSE
!                 V( (i-1)*nl+il, (j-1)*nl+jl ) =  0.0D+00
!              ENDIF
                 
                    
!                 ENDIF


           ENDDO
        ENDDO
!
     ENDDO bsp_j
  ENDDO bsp_i
!
!
  IF(print_out=='yes') THEN
     CALL print_mx(dim_h_1, v, 'v', 'e')
  ENDIF

END SUBROUTINE MAKE_V_MOLECULAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (Magnetized hydrogen) 
SUBROUTINE MAKE_V_MAGNETIC(T, V)

  USE ang,         only : cleb
  USE utils,       only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE
!
!declarations!
!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v 
!
! locals !
!
  REAL(DPK), DIMENSION(nb-1,nb-1)    :: mx_r0, mx_r2
  REAL(DPK), DIMENSION(nl,nl)        :: vlm
  INTEGER i,j,il,jl
  INTEGER dim_h_1, dim_h_2
  INTEGER la, lb
  REAL(DPK) lab
  REAL(dpk) r2, r0
  REAL(DPK) h_b_o, h_l_m, h_b_m
 ! EXTERNAL f

!executable part!....
  WRITE(*,*) '# h1e::make_v_magnetic:       call.'
!  stop

  CALL  mat_bsp( nb, t,   mx_r0,  f, 0, 0 )   ! < B_i | x^{0} | B_j  >    
  CALL  mat_bsp( nb, t,   mx_r2,  f, 2, 0 )    ! < B_i | x^(2) | B_j  >    

!
  dim_h_1 = nl* SIZE(mx_r0, 1)
  dim_h_2 = nl* SIZE(mx_r0, 1)
!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print O,B:"
     CALL print_mx(nb-1, mx_r0,  'B', 'f')
     CALL print_mx(nb-1, mx_r2,  'O', 'e')
  ENDIF


! Atoms in strong magnetic fields, Ruder, Wunner, Herold, Geyer, Springer, Berlin, 1994
! magnetic contribution to field-free hydrogenic hamiltonian 
! l's-coupling          : h_l_m
!  field-free hydrogenic states are recovered if  b_f == 0


!
!.... calculate vlm first (independent on B-splines basis)
!

  la = l0
  DO il = 1, nl
     lb = l0
     DO jl = 1, nl
        
        lab = (2.0D+00 * DBLE(la) + 1)*(2.0d+00 * DBLE(lb) + 1)
        
        vlm(il,jl) = SQRT(lab)                                      &
                   * cleb( 2*la,     0, 4, 0, 2*lb,   0   )         &
                   * cleb( 2*la, 2*m_l, 4, 0, 2*lb, 2*m_l )
        
       
        lb = lb + 2
     ENDDO
     la = la + 2
  ENDDO

!
!...
!

  V = 0.0d+00
  

  DO i = 1, nb-1
     DO j =1 , nb-1

        r0 = mx_r0(i,j)
        r2 = mx_r2(i,j)
                 
        la = l0
        DO il = 1, nl

           lb = l0

              DO jl = 1, nl

                 IF( il == jl ) THEN


                    V( (i-1)*nl+il, (j-1)*nl+jl) = r2 * b_f**2 /3.0D+00   + r0 * b_f* (m_l - 1)

              ELSE

!                 lab = (2.0D+00 * DBLE(la) + 1)*(2.0d+00 * DBLE(lb) + 1)
!                 h_l_m = sqrt(lab)                                      &
!                       * cleb( 2*la,     0, 4, 0, 2*lb,   0   )         &
!                       * cleb( 2*la, 2*m_l, 4, 0, 2*lb, 2*m_l )
!                 V( (i-1)*nl+il, (j-1)*nl+jl) =  - r2 * h_l_m * b_f**2 /3.0D+00 

                 V( (i-1)*nl+il, (j-1)*nl+jl) =  - r2 * vlm(il,jl) * b_f**2 /3.0D+00 

              ENDIF  
              lb = lb + 2
           ENDDO
           la = la + 2
        ENDDO
 
     ENDDO
  ENDDO


  IF(print_out=='yes') THEN
     CALL print_mx(dim_h_1, v, 'v', 'e')
  ENDIF

END SUBROUTINE MAKE_V_MAGNETIC
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx B
SUBROUTINE MAKE_B(T, B)

  USE utils,      only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE
!declarations!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: b 
! locals !
  REAL(DPK), DIMENSION(nb-1,nb-1)        :: mx_r0
  INTEGER i,j,il,jl
  INTEGER dim_h_1, dim_h_2
  INTEGER la, lb
  REAL(dpk) r0
!  external f
!executable part!....


  CALL  mat_bsp( nb, t,  mx_r0,  f, 0, 0 )   ! < B_i | r^0 | B_j  > = <B_i|B_j>   

  WRITE(*,*) '# h1e::make_B:       call.'
!
  dim_h_1 = nl* SIZE(mx_r0, 1)
  dim_h_2 = nl* SIZE(mx_r0, 1)
!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print B:"
     CALL print_mx(nb-1, mx_r0,  'B', 'f')
  ENDIF

!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES


! Construct the  nl*(nb-1) x nl*(nb-1) overlap matrix  BLMAX(il;jl) from  B(nb-1,nb-1)


  B = 0.0d+00
  
  DO i = 1, nb-1
     DO j =1 , nb-1

        r0 = mx_r0(i,j)
                 
        la = l0
        DO il = 1, nl

           lb = l0

              DO jl = 1, nl

                 IF( il == jl ) THEN

                 B( (i-1)*nl+il, (j-1)*nl+jl) =  r0

              ELSE

                 B( (i-1)*nl+il, (j-1)*nl+jl) =  0.0D+00  
              ENDIF  

              lb = lb + 2
           ENDDO
           la = la + 2
        ENDDO
 
     ENDDO
  ENDDO


  IF(print_out=='yes') THEN
     CALL print_mx(dim_h_1, b, 'b', 'e')
  ENDIF

END SUBROUTINE MAKE_B
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE MAKE_H0(T, H0)

  USE utils,      only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE

!declarations!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: h0
! locals !
  REAL(DPK), DIMENSION(nb-1,nb-1)        :: mx_dr_0
  REAL(DPK), DIMENSION(nb-1,nb-1)        :: mx_r_1, mx_r_2
  INTEGER i,j,il,jl
  INTEGER dim_h_1, dim_h_2
  INTEGER la, lb
  REAL(DPK) lab
  REAL(dpk) p2, r_2, r_1
!  external f

!!!....  

  CALL  mat_bsp( nb, t,  mx_dr_0,  f,  0, 1 )   ! < B'_i| 1      | B'_j >     kinetic     operator
  CALL  mat_bsp( nb, t,   mx_r_1,  f, -1, 0 )   ! < B_i | x^{-1} | B_j  >    coulombing  potential
  CALL  mat_bsp( nb, t,   mx_r_2,  f, -2, 0 )   ! < B_i | x^(-2) | B_j  >    centrifugal operator

!executable part!....

  WRITE(*,*) '# h1e::make_h:       call.'
!  stop
!
  dim_h_1 = nl* SIZE(mx_dr_0, 1)
  dim_h_2 = nl* SIZE(mx_dr_0, 1)
!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print bp,v,c:"
     CALL print_mx(nb-1, mx_dr_0, 'BP', 'e')
     CALL print_mx(nb-1, mx_r_1,   'V', 'e')
     CALL print_mx(nb-1, mx_r_2,   'C', 'e')
  ENDIF


! Construct the  nl*(nb-1) x nl*(nb-1) hamiltonian         H(il;jl) from h0(nb-1,nb-1), O,A

  WRITE(*,*) '# make_h0::         b_f = ', b_f
  WRITE(*,*) '# make_h0::          l0 = ', l0


  H0 = 0.0d+00
  
  bsp_i: DO i = 1, nb-1
     bsp_j: DO j =1 , nb-1

        p2   = mx_dr_0(i,j)
        r_1  = mx_r_1(i,j)
        r_2  = mx_r_2(i,j)

        la = l0
        DO il = 1, nl

           lb = l0

              DO jl = 1, nl

                 IF( il == jl ) THEN
                    H0( (i-1)*nl+il, (j-1)*nl+jl) = p2  /(2.0D+00*ma)         &
                                                  + r_2 * la*(la+1)/2.0D+00   &
                                                  - r_1 * za                   
                 ELSE
                    H0( (i-1)*nl+il, (j-1)*nl+jl) = 0.0D+00
              ENDIF  
              lb = lb + 2
           ENDDO
           la = la + 2
        ENDDO
        
     ENDDO bsp_j
  ENDDO bsp_i

  IF(print_out=='yes') THEN
     CALL print_mx(dim_h_1, h0, 'h0', 'e')
  ENDIF

END SUBROUTINE MAKE_H0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DIAGONALIZE_F02FHF(H,B)

  USE ang,        only : cleb
  USE utils,      only : print_mx
!......................

  IMPLICIT NONE

!declarations!
  REAL(DPK), DIMENSION(:,:), intent(in):: h         ! full matrix
  REAL(DPK), DIMENSION(:,:), intent(in):: b         ! full matrix

! locals !
  INTEGER I, J
  INTEGER N, K

!F02FHF (eigenvalues only)
  REAL(DPK), DIMENSION( NL*KB, NL*(NB-2) )                        :: hb, bb    !banded matrices
  REAL(DPK), DIMENSION( MAX(NL*(nb-2),4*(NL*kB-1)*(2*NL*KB-1) ) ) :: w
  REAL(DPK), DIMENSION( NL*(NB-2) )                               :: en
  INTEGER LWORK, INFO 

!executable part!....

  WRITE(*,*) '# h1e::solve_se_eigensystem_nag:       call.'


!...

  K   = NL *  KB - 1
  N   = NL * (NB - 2)

!...

!
! hb, bb, DIMENSION(  NL*KB, NL*(NB-2)) 
!
! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
!

  HB = 0.0D+00
  BB = 0.0D+00

  DO J = 1, N
     DO I = MAX(1, J - K), J

        HB( K + 1 + I - J, J ) =   H(I, J)
        BB( K + 1 + I - J, J ) =   B(I, J)
     ENDDO
  ENDDO


  !  diagonalize    :
  !  FO2FHF/NAG/M20 :     A * C = e * B * C,   where, A,B banded real symmetric
  !


  INFO = 0
  LWORK = MAX( NL*(nb-2), 4*(NL*kB-1)*(2*NL*KB-1))    

!...

  CALL F02FHF(N, K, HB, K+1, K, BB, K+1, EN, W, LWORK, INFO)

!...
      
  IF(INFO.NE.0) THEN 
     
     WRITE(*,*) '# ERROR IN  F02FHF/NAG/M20 SUBROUTINE,    INFO = ', INFO 
     WRITE(*,*) '#'
     STOP
     
  ELSE  

     OPEN(10, file='out/en-nag.out') 
     DO  I = 1 , NL*(NB-2) 
        
        WRITE(10,"(I4,3X,3E20.10)") I, EN(I)
        
     ENDDO
     CLOSE(10)
  ENDIF
      
END SUBROUTINE DIAGONALIZE_F02FHF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DIAGONALIZE_DSBGV(B, H0, V, EN, CE)

  !declarations!
  IMPLICIT NONE
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: H0
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: B
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: V
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT) :: CE
  REAL(DPK), DIMENSION(:),    INTENT(INOUT) :: EN

  !locals!
  INTEGER I, J
  INTEGER N, K

  ! DSBGV/LAPACK
  REAL(DPK), DIMENSION(  NL*KB, NL*(NB-2) )     :: HB, BB
  REAL(DPK), DIMENSION(  NL*(NB-2), NL*(NB-2) ) :: Z      !dsbgv
  REAL(DPK), DIMENSION( 3*NL*(NB-2) )           :: w      !dsbgv
  INTEGER INFO 

!executable part!

!...

  WRITE(*,*) '# h1e::diagonalize_dsbgv:   call.'

!...

  K   = NL *  KB - 1
  N   = NL * (NB - 2)

!............................

! hb, bb, DIMENSION(  NL*KB, NL*(NB-2)) 
!
! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
!

      HB = 0.0D+00
      BB = 0.0D+00

      DO J = 1, N
         DO I = MAX(1, J - K), J
            HB( K + 1 + I - J, J ) =   H0(I, J) + V(I,J)
            BB( K + 1 + I - J, J ) =   B(I, J) 
         ENDDO
      ENDDO

      !
      !  DSBGV/LAPACK :     A * C = e * B * C,    where, A,B banded real symmetric
      !
      !   EV NORMALIZED AS : Z^T * C * Z = 1
      !   EN(I) <====> Z(J,I)
      !

      INFO = 0

      CALL DSBGV('V', 'U', N, K, K, HB, K+1, BB, K+1, EN, Z, N, W, INFO)
      
      !...
      IF(INFO.NE.0) THEN 

         WRITE(*,*) '# ERROR IN  DSBGV/LAPACK SUBROUTINE,    INFO = ', INFO 
         WRITE(*,*) '#'
         STOP
      ELSE  

         OPEN(10, file='out/en-dsbgv.out') 
         DO  I = 1 , NL*(NB-2) 

            WRITE(10,"(I4,3X,3E20.10)") I, EN(I)

         ENDDO
         CLOSE(1)

         CE = TRANSPOSE(Z)                 !assign Z(NB-2,NB-2) ---> CE(NB-2,NB-1)
         CE(:, NL*(NB-1))  = 0.0D+00

         DO i = 1, nl*(nb-2)
            WRITE(*,*) i, -en(i), -2*en(i)
         ENDDO
      ENDIF
     
      
    END SUBROUTINE DIAGONALIZE_DSBGV
!f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90f90
    SUBROUTINE LINEARIZE_F07BEF(H, B, EN, CE)  
!  USE nag_f77_f_chapter
  USE utils, ONLY : check
!  USE utils, ONLY : print_mx, check

!declarations!
  IMPLICIT NONE
  REAL(DPK), DIMENSION(:,:),           INTENT(IN) :: H
  REAL(DPK), DIMENSION(:,:),           INTENT(IN) :: B
  REAL(DPK), DIMENSION(:),          INTENT(inout) :: en
  REAL(DPK), DIMENSION(:,:),        INTENT(inout) :: ce

! locals !
  INTEGER i, j
!F07BEF!
  INTEGER K, N
  REAL(DPK), DIMENSION( 3*NL*KB-2, NL*(NB-1)) :: A_E  ! linear
  REAL(DPK), DIMENSION( NL*(NB-1),    1     ) :: C0
  INTEGER,   DIMENSION(   NL*(NB-1) )         :: IW
  INTEGER INFO
!
  REAL(dpk) norm
  INTEGER l, ie
  INTEGER nout
  INTEGER dim, dim_1, dim_2  
  integer er

!executable!
!...
  WRITE(*,*) '# h1e::linearize_f07bef:                  call.'

!...

  NOUT     = 12

!...

  K  = NL *  KB - 1          ! number of sub-diagonials (upper diagonals)
  N  = NL * (NB - 1)         ! dimension of coefficient vector

!...

  er = 1


!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES


  OPEN(NOUT, FILE='out/subsolve_lin.out')

  WRITE(nout, *) '#       LINEAR SOLUTION : '

!...


  dim = size(en)
  dim_1 = SIZE(ce,1)
  dim_2 = SIZE(h, 1)

  WRITE(*,*) '# nof states          ns = ', dim_1, nbs + ncs 
  WRITE(*,*) '# nof basis           nc = ', dim_2, nl * ( nb - 1 )

  CALL check(dim_1, dim, er)
  CALL check(dim_2, nl*(nb-1), er+1)

    
  DO ie = 1, dim                      ! solve for 'ns' states: en(ie), ie = 1, ..., ns

!...           A(E) = H - E*B 
 
     DO J = 1, N
        DO I = MAX(1, J - K), MIN(N, J + K)
           
           A_E(2*K + I - J + 1, J) = H(I,J) - EN(IE) * B(I,J)         

        ENDDO
     ENDDO


!... solve     A(E) * C = C0,

     IW = 0.0D+00

     CALL F07BDF(N, N, K, K, A_E, 3*NL*KB-2, IW, INFO)      ! factorize

     IF(INFO.NE.0) THEN 

        WRITE(*,*) '# subsolve:: f07bdf: Error in factorization routine.'
        WRITE(*,*) '# subsolve:: f07bdf:                         info = ', info

        CALL X04CEF(N, N, K, 2*K, A_E, 3*NL*KB-2,'X04CEF: A=PLU, U:', INFO)

     ELSE

        C0 = 0.0D+00                        ! initialize C0
        C0(NB-1, 1) = 1.0D+00

        IF(ie.EQ.1) THEN
           WRITE(*,*) '# c0 = ', ie,  SIZE(c0, 1), SIZE(c0, 2) 
        ENDIF


        CALL F07BEF('N', N, K, K, 1, A_E, 3*NL*KB-2, IW, C0, N, INFO)

          ENDIF


          IF(INFO.NE.0) THEN 
             WRITE(nout,*) '#         ERROR IN   F07BEF          '
             WRITE(nout,*) '#'
             WRITE(nout,*) '#                     INFO = ', INFO
             WRITE(nout,*) '#'
             STOP
          ELSE

             WRITE(nout, *) '#   INPUT     ENERGY      E = ', EN(IE)
             WRITE(nout, *) '#'
             
                          
             NORM = DOT_PRODUCT( C0(:,1), MATMUL( B, C0(:, 1) ))


!             CE(IE,:) = C(:,1) /SQRT(NORM) 
             

             ! transform output vector from c_i(l=0,1,..,nl) ---> c_l(i=1,2,..,nb-1)

             DO l = 1, nl
                DO j = 1, nb-1 

                   CE(ie,j +(l-1)*(nb-1)) = C0(l +(j-1)*nl,1) /DSQRT(NORM) 

                ENDDO
             ENDDO

          ENDIF

          WRITE(NOUT, *) IE, EN(IE) 

          IF(IE.EQ.30.AND.L.EQ.0) THEN
             DO  J = 1, (NB-1)/2
               
                WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, CE(IE,J), J+NB/2, CE(IE, J+NB/2)
             ENDDO
          ENDIF

       ENDDO

       CLOSE(NOUT)
       
     END SUBROUTINE LINEARIZE_F07BEF
    !###########################################
    SUBROUTINE ENERGY_SPECTRUM(L, EN, MODE)
 
      USE param
      USE DATA, ONLY: write_v
      
      IMPLICIT NONE

      INTEGER L
      REAL(DPK), DIMENSION(:), INTENT(INOUT) :: en
      INTEGER NS
      CHARACTER(LEN=*)  MODE
!... locals
      INTEGER I

!.....................

      WRITE(*,*) '# energy_spectrum:    size of en = ', size(en)

      ns = SIZE(en)
      
      IF (MODE =='OUT') THEN 

         DO  I = L + 1 , NS
            
            IF(I.LE.NBS + L) THEN 

               EN(I) = - 0.5D+00 * (za/I)**2

            ELSE IF(I.LE.NS) THEN 
               
               EN(I) =   EN_1E  +    DE_1E * ( I - NBS - 1 - L) 
     
            ENDIF

         ENDDO

         CALL write_v(l, en, "en1e-")

      ELSE IF (MODE =='INOUT') THEN      ! do not store

         DO  I = 1 , NS
            EN(I) =   EN_1E  +    DE_1E * ( I - 1) 
         ENDDO

      ELSE IF (MODE == 'IN') THEN 

         CALL write_v(l, en, "en1e-")
      ENDIF
    END SUBROUTINE ENERGY_SPECTRUM


  END PROGRAM H1EZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF
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



